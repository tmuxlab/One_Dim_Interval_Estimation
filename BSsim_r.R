# 必要パッケージ（未インストールなら一度だけ）
pkgs <- c("peakRAM", "readr", "tibble", "fs")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 共通：ログ出力ヘルパ ----
.ensure_dir <- function(path) { if (!fs::dir_exists(path)) fs::dir_create(path) }

.write_log_row <- function(log_path, row_df) {
  header_needed <- !file.exists(log_path)
  readr::write_csv(row_df, log_path, append = !header_needed)
}

run_and_log <- function(fun, fun_name,
                        Estimate, roopj,
                        EstimationMethod, thresholdRule, thresholdValue,
                        directory_path, DS,
                        log_path) {
  t0 <- proc.time()
  
  fit_ok <- TRUE
  res <- NULL
  pr <- peakRAM::peakRAM({
    res <- tryCatch(
      fun(Estimate, roopj, EstimationMethod, thresholdRule, thresholdValue, directory_path),
      error = function(e) {
        fit_ok <<- FALSE
        NULL
      }
    )
  })
  
  # fun() からの返り値で fit_ok / params を拾う
  if (!is.null(res) && is.list(res)) {
    if (!is.null(res$fit_ok)) fit_ok <- isTRUE(res$fit_ok)
  }
  # パラメータ列の標準化
  param_names <- c("a","b","mu0","sita","beta")
  param_vals  <- setNames(as.list(rep(NA_real_, length(param_names))), param_names)
  if (!is.null(res) && is.list(res) && !is.null(res$params)) {
    got <- intersect(names(res$params), param_names)
    for (nm in got) param_vals[[nm]] <- suppressWarnings(as.numeric(res$params[[nm]]))
  }
  
  t1 <- proc.time() - t0
  wall_sec <- suppressWarnings(as.numeric(pr$Elapsed_Time_sec[1]))
  if (is.na(wall_sec)) wall_sec <- as.numeric(t1[["elapsed"]])
  peak_mb <- suppressWarnings(as.numeric(pr$Peak_RAM_Used_MiB[1]))
  
  row <- tibble::tibble(
    timestamp    = as.character(Sys.time()),
    dataset      = DS,
    J            = roopj,
    transform    = EstimationMethod,
    threshold    = thresholdRule,
    mode         = thresholdValue,
    user_sec     = as.numeric(t1[["user.self"]]),
    system_sec   = as.numeric(t1[["sys.self"]]),
    elapsed_sec  = wall_sec,
    peak_ram_mib = peak_mb,
    method       = fun_name,
    pid          = Sys.getpid(),
    fit_ok       = fit_ok,
    param_a      = param_vals$a,
    param_b      = param_vals$b,
    param_mu0    = param_vals$mu0,
    param_sita   = param_vals$sita,
    param_beta   = param_vals$beta
  )
  
  .write_log_row(log_path, row)
  invisible(row)
}

# ---- メイン：BSSim（計測版） ----
BSSim <- function(DS, roopJ, directory_path, out_path,
                  EstimationMethod, thresholdRule, thresholdValue) {
  
  base_dir <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), directory_path) # BS-output 側
  .ensure_dir(base_dir)
  
  # ログは BS-output 側へ
  LOG_PATH <- file.path(base_dir, "method_bench_log.csv")
  
  # ★ 追加：実行する手法（必要に応じてON/OFF）
  methods <- list(
    #COS  = COSsimulate,
    #GSm  = GSmsimulate,
    #GSp  = GSpsimulate,
    #Ord  = Ordsimulate,
    #Thi  = Thisimulate,
    # RGS = RGSsimulate,
    # Tim = Timsimulate,
    # Int = Intsimulate,
    ThiInt = ThiIntsimulate
  )
  
  for (roopj in 2:roopJ) {
    # 推定結果は out_path（前段の出力＝今回の入力）から読む
    dataPath <- paste0(
      dirname(rstudioapi::getSourceEditorContext()$path),
      out_path, "k", EstimationMethod, "-", thresholdRule, "-", thresholdValue,
      "(J=", roopj, ").txt"
    )
    
    # ★ 推定結果ファイルの存在チェック（落ちやすいので丁寧に）
    if (!file.exists(dataPath)) {
      warning(sprintf("Estimate file not found: %s", dataPath))
      next
    }
    
    row_tbl <- read.table(dataPath)
    estimate <- row_tbl[[1]]
    Estimate <- as.numeric(estimate[2:nrow(row_tbl)])
    
    # ★ 追加：各手法を計測してログに追記
    for (nm in names(methods)) {
      fun <- methods[[nm]]
      try({
        run_and_log(fun, nm, Estimate, roopj,
                    EstimationMethod, thresholdRule, thresholdValue,
                    directory_path = directory_path, DS = DS,
                    log_path = LOG_PATH)
      }, silent = TRUE)
    }
  }
}
