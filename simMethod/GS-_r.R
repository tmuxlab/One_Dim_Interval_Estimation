GSmsimulate = function(estimate, roopj, estimationMethod, thresholdRule, thresholdValue, directory_path)
{
  tn <- length(estimate)
  timeindex <- 1:tn
  x <- c(0, timeindex)
  y <- c(0, estimate)
  
  row <- 10000
  col <- 500
  arrivals <- matrix(0, nrow = row, ncol = col)
  counts   <- matrix(0, nrow = row, ncol = tn)
  
  # --- 推定（nls）を安全化：失敗時は log-linear にフォールバック ---
  fit_ok <- TRUE
  GSm <- tryCatch(
    {
      nls(y ~ a * exp(b * x),
          start = list(a = 1, b = -10),
          algorithm = "port",
          lower = list(a = 0.0005, b = -100),
          upper = list(a = 1000,   b = -0.0005))
    },
    error = function(e) {
      fit_ok <<- FALSE
      y_pos <- pmax(y, 1e-8)
      co <- coef(lm(log(y_pos) ~ x))
      list(a = exp(co[[1]]), b = min(co[[2]], -0.0005))
    }
  )
  
  if (isTRUE(fit_ok)) {
    a_hat <- as.numeric(coef(GSm)[["a"]])
    b_hat <- as.numeric(coef(GSm)[["b"]])
  } else {
    a_hat <- GSm$a
    b_hat <- GSm$b
  }
  # 健全化
  if (!is.finite(a_hat) || a_hat <= 0) a_hat <- 1e-6
  if (!is.finite(b_hat)) b_hat <- -1e-6
  
  sita <- -a_hat / b_hat
  beta <- -b_hat
  if (!is.finite(sita) || sita < 0) { sita <- max(sita, 0); fit_ok <- FALSE }
  if (!is.finite(beta) || beta <= 0) { beta <- 1e-6; fit_ok <- FALSE }
  
  get_nhpp_realization <- function(sita, beta, tn){
    i <- 1; t <- 0
    m <- rpois(1, sita)
    if (m <= 0) return(numeric(0))
    A <- numeric(0)
    while (i <= m) {
      t <- rexp(1, rate = 1) / beta / (m - i + 1) + t
      if (t > tn) break
      A <- c(A, t)
      i <- i + 1
    }
    A
  }
  
  brks <- 0:tn
  for (i in 1:row){
    t_samples <- get_nhpp_realization(sita, beta, tn)
    if (length(t_samples) && t_samples[length(t_samples)] > tn) {
      t_samples <- t_samples[t_samples <= tn]
    }
    ln <- length(t_samples); ln <- min(ln, col)
    arrivals[i, ] <- NA
    if (ln > 0) arrivals[i, 1:ln] <- t_samples[1:ln]
    counts[i, ] <- hist(t_samples, plot = FALSE, breaks = brks, right = FALSE)$counts
  }
  
  base <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), directory_path)
  write.table(arrivals, paste0(base, "BS-GS--", estimationMethod, "-", thresholdRule, "-", thresholdValue, "(J=", roopj, ").dat"),
              append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(counts,   paste0(base, "BSCount-GS--", estimationMethod, "-", thresholdRule, "-", thresholdValue, "(J=", roopj, ").dat"),
              append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  return(invisible(list(
    fit_ok = fit_ok,
    params = list(
      a    = a_hat,                  # nls 推定 a
      b    = b_hat,                  # nls 推定 b
      sita = sita,                   # 変換後
      beta = beta                    # 変換後
      # mu0 はこの手法では未定義 → 省略（NA）
    )
  )))
}
