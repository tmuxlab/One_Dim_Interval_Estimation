# ============================================================
# 共通ユーティリティ：指数モデル y ≈ a * exp(b * x) の初期値をデータから推定
# ------------------------------------------------------------
# - y<=0 は eps でクリップして log を安定化
# - log-linear 回帰の傾きをベースに b を作成（必ず負にクリップ）
# - a は exp(切片) をベースに、非有限/非正を防ぐ
# ============================================================
get_init_ab <- function(x, y, eps = 1e-8) {
  y_clip <- pmax(y, eps)
  fit <- try(suppressWarnings(lm(log(y_clip) ~ x)), silent = TRUE)
  
  if (inherits(fit, "try-error")) {
    a0 <- 1
    b0 <- -1e-3
  } else {
    co <- coef(fit)
    a0 <- exp(as.numeric(co[[1]]))
    b0 <- as.numeric(co[[2]])
    if (!is.finite(a0) || a0 <= 0) a0 <- 1
    if (!is.finite(b0)) b0 <- -1e-3
  }
  # 物理制約：減衰（b<0）を必ず満たすように弱くクリップ
  if (b0 >= -1e-8) b0 <- -1e-3
  c(a = a0, b = b0)
}

# ============================================================
# COSsimulate（改良版）— 出力形式は Ordsimulate と同様
#   arrivals: m × Lmax（不足は NA パディング）
#   counts  : m × t_n
# ============================================================
COSsimulate = function(estimate, roopj, estimationMethod, thresholdRule, thresholdValue, directory_path)
{
  t_n <- length(estimate)
  
  # x, y 準備（従来どおり 0 を含める）
  timeindex <- 1:t_n
  x <- c(0, timeindex)
  y <- c(0, estimate)
  
  # --- 初期値をデータから自動推定 ---
  init <- get_init_ab(x, y)
  
  # --- nls を制約付きで実行（失敗時はフォールバック） ---
  fit_ok <- TRUE
  COS <- tryCatch(
    {
      nls(
        y ~ a * exp(b * x),
        start = list(a = init["a"], b = init["b"]),
        algorithm = "port",
        lower = list(a = 0.000005, b = -100),
        upper = list(a = 1000,     b = -0.00005)
      )
    },
    error = function(e) {
      fit_ok <<- FALSE
      y_pos <- pmax(y, 1e-8)
      lmfit <- lm(log(y_pos) ~ x)
      co <- coef(lmfit)
      list(a = exp(co[[1]]), b = min(co[[2]], -0.00005))
    }
  )
  
  if (fit_ok) {
    lambda <- as.numeric(coef(COS)[["a"]])
    alpha1 <- as.numeric(coef(COS)[["b"]])
  } else {
    lambda <- COS$a
    alpha1 <- COS$b
  }
  
  # --- パラメータ健全化 ---
  if (!is.finite(alpha1) || alpha1 >= -1e-8) alpha1 <- -1e-8
  if (!is.finite(lambda) || lambda <= 0)     lambda <- 1e-6
  
  # 累積強度の総量（Poisson 平均）：指数モデルの厳密式
  C    <- exp(alpha1 * t_n) - 1
  myu0 <- lambda * C / alpha1  # = Λ(t_n)
  
  # 総強度が非正/非有限なら全ゼロで返す
  if (!is.finite(myu0) || myu0 <= 0) {
    m <- 10000
    counts   <- matrix(0,  nrow = m, ncol = t_n)
    arrivals <- matrix(NA, nrow = m, ncol = 1)
    base <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), directory_path)
    write.table(arrivals, paste0(base, "BS-COS-", estimationMethod, "-", thresholdRule, "-", thresholdValue,
                                 "(J=", roopj, ").dat"),
                append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE)
    write.table(counts,   paste0(base, "BSCount-COS-", estimationMethod, "-", thresholdRule, "-", thresholdValue,
                                 "(J=", roopj, ").dat"),
                append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE)
    return(invisible(NULL))
  }
  
  # --- NHPP(指数) 逆変換で生成（従来ロジック／数値安定だけ強化） ---
  get_nhpp_realization <- function(myu0, alpha1, C) {
    n <- rpois(1, myu0)
    if (n <= 0) return(numeric(0))
    U <- sort(runif(n))
    T <- log1p(C * U) / alpha1   # log(1+z) で安定化
    T
  }
  
  # --- 出力形式は Ordsimulate と同様 ---
  m <- 10000
  arrivals_list <- vector("list", m)
  counts <- matrix(0, nrow = m, ncol = t_n)
  brks <- 0:t_n
  
  for (i in 1:m) {
    t_samples <- get_nhpp_realization(myu0, alpha1, C)
    
    # 念のための範囲ガード（理論上 tn 超えは出ない）
    if (length(t_samples) && t_samples[length(t_samples)] > t_n) {
      t_samples <- t_samples[t_samples <= t_n]
    }
    
    arrivals_list[[i]] <- t_samples
    counts[i, ] <- hist(t_samples, plot = FALSE, breaks = brks)$counts
  }
  
  # 可変長 → NA パディング
  max_length <- max(1, max(sapply(arrivals_list, length)))
  arrivals <- matrix(NA, nrow = m, ncol = max_length)
  for (i in 1:m) {
    li <- length(arrivals_list[[i]])
    if (li > 0) arrivals[i, 1:li] <- arrivals_list[[i]]
  }
  
  # --- 出力 ---
  base <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), directory_path)
  write.table(arrivals, paste0(base, "BS-COS-", estimationMethod, "-", thresholdRule, "-", thresholdValue,
                               "(J=", roopj, ").dat"),
              append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(counts,   paste0(base, "BSCount-COS-", estimationMethod, "-", thresholdRule, "-", thresholdValue,
                               "(J=", roopj, ").dat"),
              append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE)

  return(invisible(list(
    fit_ok = fit_ok,                 # 既に持っているフラグ
    params = list(
      a   = lambda,                  # nls 推定 a
      b   = alpha1,                  # nls 推定 b
      mu0 = myu0                     # Λ(t_n)
      # sita, beta は COS では未使用なので省略（NA で入ります）
    )
  )))
}
