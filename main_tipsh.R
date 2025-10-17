# ================== Tipsh 用 メトリクス記録ユーティリティ ==================
# できるだけ正確に現在メモリ(MB)を取得
get_mem_mb_tipsh <- function() {
  if (requireNamespace("lobstr", quietly = TRUE)) {
    return(as.numeric(lobstr::mem_used()) / 1024^2)
  }
  if (.Platform$OS.type == "windows" && exists("memory.size")) {
    return(as.numeric(memory.size()))
  }
  if (file.exists("/proc/self/status")) {
    vmrss <- tryCatch(system("grep VmRSS /proc/self/status | awk '{print $2}'", intern = TRUE),
                      error = function(e) NA)
    if (!is.na(vmrss) && length(vmrss) > 0) return(as.numeric(vmrss[1]) / 1024) # kB→MB
  }
  return(NA_real_)
}

# Tipsh 実行を計測してCSVに1行追記
run_tipsh_logged <- function(ds, threshold_mode, var, J_index, dataset_id, log_path) {
  mem_before <- get_mem_mb_tipsh()
  tm <- system.time({
    res <- Tipsh(ds, threshold_mode, var, J_index)  # ← そのままTipshを呼ぶ
  })
  mem_after <- get_mem_mb_tipsh()
  
  log_row <- data.frame(
    timestamp       = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    dataset         = dataset_id,
    J               = J_index,
    transform       = "TI",          # Tipsh = translation-invariant
    threshold       = "ldt",         # Tipsh内部で固定
    mode            = threshold_mode, # "h" or "s"
    user_sec        = unname(tm[["user.self"]]),
    system_sec      = unname(tm[["sys.self"]]),
    elapsed_sec     = unname(tm[["elapsed"]]),
    mem_before_mb   = mem_before,
    mem_after_mb    = mem_after,
    mem_delta_mb    = if (is.na(mem_before) || is.na(mem_after)) NA_real_ else (mem_after - mem_before),
    stringsAsFactors = FALSE
  )
  
  dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)
  write.table(
    log_row, file = log_path, sep = ",", row.names = FALSE,
    col.names = !file.exists(log_path), append = TRUE, quote = TRUE
  )
  return(res)
}
# ===========================================================================

# Tipsh用ログの保存先（Wseとは別ファイルに）
LOG_PATH_TIPSH <- "./output/metrics/tipsh_metrics.csv"

# ---- ここから Tipsh セクション（コメントアウトを外して置き換え）----
# TIPSH: Translation-Invariant（データ変換なし、しきい値は ldt 固定）
directory_path = "./output/NDT_WSE/TI/"
dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)

for (j in 1:4) {
  dataPath = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/DS/DS",j,".txt")
  ds_tbl = read.table(dataPath)[2]
  ds = as.numeric(ds_tbl$V2)
  
  for (i in 2:log(GetGroupLength(length(ds)), base = 2)) {
    # hard / soft をそれぞれ実行＆ログ追記
    hard = run_tipsh_logged(ds, "h", 1, i, j, LOG_PATH_TIPSH)
    soft = run_tipsh_logged(ds, "s", 1, i, j, LOG_PATH_TIPSH)
    
    # 出力（従来のファイル書き出しと同等）
    edata = list(
      hard = round(hard$EstimationData, digits = 3),
      soft = round(soft$EstimationData, digits = 3)
    )
    txt_path = create_txt(i, j, directory_path, "ldt", "TI")
    write.table(edata$hard, txt_path$hdata, row.names = FALSE)
    write.table(edata$soft, txt_path$sdata, row.names = FALSE)
  }
}
# ---- ここまで Tipsh セクション ----
