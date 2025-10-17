# install packages "tidyverse"
# install.packages("tidyverse")
suppressPackageStartupMessages(library(tidyverse))

# ========= 追加: 計測ユーティリティ =========
# メモリ(MB)をなるべく正確に取得
get_mem_mb <- function() {
  if (requireNamespace("lobstr", quietly = TRUE)) {
    return(as.numeric(lobstr::mem_used()) / 1024^2)
  }
  if (.Platform$OS.type == "windows" && exists("memory.size")) {
    return(as.numeric(memory.size()))  # MB
  }
  if (file.exists("/proc/self/status")) {
    vmrss <- tryCatch(system("grep VmRSS /proc/self/status | awk '{print $2}'", intern = TRUE),
                      error = function(e) NA)
    if (!is.na(vmrss) && length(vmrss) > 0) return(as.numeric(vmrss[1]) / 1024)  # kB->MB
  }
  return(NA_real_)
}

# Wse を包んで計測＆CSV追記するラッパ
run_wse_logged <- function(ds, transform, threshold_name, threshold_mode, var, J_index,
                           dataset_id, log_path) {
  mem_before <- get_mem_mb()
  tm <- system.time({
    res <- Wse(ds, transform, threshold_name, threshold_mode, var, J_index)
  })
  mem_after <- get_mem_mb()
  
  log_row <- data.frame(
    timestamp       = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    dataset         = dataset_id,
    J               = J_index,
    transform       = transform,
    threshold       = threshold_name,
    mode            = threshold_mode,
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
# ===========================================

# ログ出力先
LOG_PATH <- "./output/metrics/wse_metrics.csv"

# Load Hal wavelet estimation module
WSE_Path = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/src/WaveletShrinkageEstimation.R")
print("Load Hal wavelet estimation module")
source(WSE_Path)

WT_Path = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/src/WaveletTransform.R")
print("Load Hal wavelet Transformation module")
source(WT_Path)

for(j in 1:4){
  # Load data set
  dataPath = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/DS/DS",j,".txt")
  ds = read.table(dataPath)[2]
  ds = as.numeric(ds$V2)
  
  # ===============================================================
  # H  (コメントアウトのままなので変更なし)
  # ===============================================================
   directory_path = "./output/NDT_WSE/"
  for(i in 2:log(GetGroupLength(length(ds)), base = 2)){
     hard = run_wse_logged(ds, "none", "ldt", "h", 1, i, j, LOG_PATH)
     soft = run_wse_logged(ds, "none", "ldt", "s", 1, i, j, LOG_PATH)
     edata = list(hard = round(hard$EstimationData, digits = 3), soft = round(soft$EstimationData, digits = 3))
     txt_path = create_txt(i,j,directory_path,"ldt","H")
     write.table(edata$hard,txt_path$hdata, row.names = FALSE)
     write.table(edata$soft,txt_path$sdata, row.names = FALSE)
  }

  # ===============================================================
  # TIPSH (コメントアウトのままなので変更なし)
  # ===============================================================
  # directory_path = "./output/NDT_WSE/TI/"
  # for(i in 2:log(GetGroupLength(length(ds)), base = 2)){
  #  hard = Tipsh(ds, "h", 1, i)
  #  soft = Tipsh(ds, "s", 1, i)
  #  edata = list(hard = round(hard$EstimationData, digits = 3), soft = round(soft$EstimationData, digits = 3))
  #  txt_path = create_txt(i,j,directory_path,"ldt","TI")
  #  write.table(edata$hard,txt_path$hdata, row.names = FALSE)
  #  write.table(edata$soft,txt_path$sdata, row.names = FALSE)
  # }
  
  # ===============================================================
  # HAT (Anscombe)
  # ===============================================================
  
  # HAT_A1
  directory_path = "./output/DT_Ans_WSE/A1/"
  dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
  for(i in 2:log(GetGroupLength(length(ds)), base = 2)){
    ut_hard = run_wse_logged(ds, "A1", "ut", "h", 1, i, j, LOG_PATH)
    ut_soft = run_wse_logged(ds, "A1", "ut", "s", 1, i, j, LOG_PATH)
    ut_edata = list(hard = round(ut_hard$EstimationData, digits = 3), soft = round(ut_soft$EstimationData, digits = 3))
    ut_txt_path = create_txt(i,j,directory_path,"ut","A1")
    write.table(ut_edata$hard,ut_txt_path$hdata, row.names = FALSE)
    write.table(ut_edata$soft,ut_txt_path$sdata, row.names = FALSE)
    
    lut_hard = run_wse_logged(ds, "A1", "lut", "h", 1, i, j, LOG_PATH)
    lut_soft = run_wse_logged(ds, "A1", "lut", "s", 1, i, j, LOG_PATH)
    lut_edata = list(hard = round(lut_hard$EstimationData, digits = 3), soft = round(lut_soft$EstimationData, digits = 3))
    lut_txt_path = create_txt(i,j,directory_path,"lut","A1")
    write.table(lut_edata$hard,lut_txt_path$hdata, row.names = FALSE)
    write.table(lut_edata$soft,lut_txt_path$sdata, row.names = FALSE)
    
    lht_hard = run_wse_logged(ds, "A1", "lht", "h", 1, i, j, LOG_PATH)
    lht_soft = run_wse_logged(ds, "A1", "lht", "s", 1, i, j, LOG_PATH)
    lht_edata = list(hard = round(lht_hard$EstimationData, digits = 3), soft = round(lht_soft$EstimationData, digits = 3))
    lht_txt_path = create_txt(i,j,directory_path,"lht","A1")
    write.table(lht_edata$hard,lht_txt_path$hdata, row.names = FALSE)
    write.table(lht_edata$soft,lht_txt_path$sdata, row.names = FALSE)
  }
  
  # HAT_A2
  directory_path = "./output/DT_Ans_WSE/A2/"
  dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
  for(i in 2:log(GetGroupLength(length(ds)), base = 2)){
    ut_hard = run_wse_logged(ds, "A2", "ut", "h", 1, i, j, LOG_PATH)
    ut_soft = run_wse_logged(ds, "A2", "ut", "s", 1, i, j, LOG_PATH)
    ut_edata = list(hard = round(ut_hard$EstimationData, digits = 3), soft = round(ut_soft$EstimationData, digits = 3))
    ut_txt_path = create_txt(i,j,directory_path,"ut","A2")
    write.table(ut_edata$hard,ut_txt_path$hdata, row.names = FALSE)
    write.table(ut_edata$soft,ut_txt_path$sdata, row.names = FALSE)
    
    lut_hard = run_wse_logged(ds, "A2", "lut", "h", 1, i, j, LOG_PATH)
    lut_soft = run_wse_logged(ds, "A2", "lut", "s", 1, i, j, LOG_PATH)
    lut_edata = list(hard = round(lut_hard$EstimationData, digits = 3), soft = round(lut_soft$EstimationData, digits = 3))
    lut_txt_path = create_txt(i,j,directory_path,"lut","A2")
    write.table(lut_edata$hard,lut_txt_path$hdata, row.names = FALSE)
    write.table(lut_edata$soft,lut_txt_path$sdata, row.names = FALSE)
    
    lht_hard = run_wse_logged(ds, "A2", "lht", "h", 1, i, j, LOG_PATH)
    lht_soft = run_wse_logged(ds, "A2", "lht", "s", 1, i, j, LOG_PATH)
    lht_edata = list(hard = round(lht_hard$EstimationData, digits = 3), soft = round(lht_soft$EstimationData, digits = 3))
    lht_txt_path = create_txt(i,j,directory_path,"lht","A2")
    write.table(lht_edata$hard,lht_txt_path$hdata, row.names = FALSE)
    write.table(lht_edata$soft,lht_txt_path$sdata, row.names = FALSE)
  }
  
  # HAT_A3
  directory_path = "./output/DT_Ans_WSE/A3/"
  dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
  for(i in 2:log(GetGroupLength(length(ds)), base = 2)){
    ut_hard = run_wse_logged(ds, "A3", "ut", "h", 1, i, j, LOG_PATH)
    ut_soft = run_wse_logged(ds, "A3", "ut", "s", 1, i, j, LOG_PATH)
    ut_edata = list(hard = round(ut_hard$EstimationData, digits = 3), soft = round(ut_soft$EstimationData, digits = 3))
    ut_txt_path = create_txt(i,j,directory_path,"ut","A3")
    write.table(ut_edata$hard,ut_txt_path$hdata, row.names = FALSE)
    write.table(ut_edata$soft,ut_txt_path$sdata, row.names = FALSE)
    
    lut_hard = run_wse_logged(ds, "A3", "lut", "h", 1, i, j, LOG_PATH)
    lut_soft = run_wse_logged(ds, "A3", "lut", "s", 1, i, j, LOG_PATH)
    lut_edata = list(hard = round(lut_hard$EstimationData, digits = 3), soft = round(lut_soft$EstimationData, digits = 3))
    lut_txt_path = create_txt(i,j,directory_path,"lut","A3")
    write.table(lut_edata$hard,lut_txt_path$hdata, row.names = FALSE)
    write.table(lut_edata$soft,lut_txt_path$sdata, row.names = FALSE)
    
    lht_hard = run_wse_logged(ds, "A3", "lht", "h", 1, i, j, LOG_PATH)
    lht_soft = run_wse_logged(ds, "A3", "lht", "s", 1, i, j, LOG_PATH)
    lht_edata = list(hard = round(lht_hard$EstimationData, digits = 3), soft = round(lht_soft$EstimationData, digits = 3))
    lht_txt_path = create_txt(i,j,directory_path,"lht","A3")
    write.table(lht_edata$hard,lht_txt_path$hdata, row.names = FALSE)
    write.table(lht_edata$soft,lht_txt_path$sdata, row.names = FALSE)
  }
  
  # ===============================================================
  # HBT (Bartlett)
  # ===============================================================
  
  # HBT_B1
  directory_path = "./output/DT_Bar_WSE/B1/"
  dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
  for(i in 2:log(GetGroupLength(length(ds)), base = 2)){
    ut_hard = run_wse_logged(ds, "B1", "ut", "h", 1, i, j, LOG_PATH)
    ut_soft = run_wse_logged(ds, "B1", "ut", "s", 1, i, j, LOG_PATH)
    ut_edata = list(hard = round(ut_hard$EstimationData, digits = 3), soft = round(ut_soft$EstimationData, digits = 3))
    ut_txt_path = create_txt(i,j,directory_path,"ut","B1")
    write.table(ut_edata$hard,ut_txt_path$hdata, row.names = FALSE)
    write.table(ut_edata$soft,ut_txt_path$sdata, row.names = FALSE)
    
    lut_hard = run_wse_logged(ds, "B1", "lut", "h", 1, i, j, LOG_PATH)
    lut_soft = run_wse_logged(ds, "B1", "lut", "s", 1, i, j, LOG_PATH)
    lut_edata = list(hard = round(lut_hard$EstimationData, digits = 3), soft = round(lut_soft$EstimationData, digits = 3))
    lut_txt_path = create_txt(i,j,directory_path,"lut","B1")
    write.table(lut_edata$hard,lut_txt_path$hdata, row.names = FALSE)
    write.table(lut_edata$soft,lut_txt_path$sdata, row.names = FALSE)
    
    lht_hard = run_wse_logged(ds, "B1", "lht", "h", 1, i, j, LOG_PATH)
    lht_soft = run_wse_logged(ds, "B1", "lht", "s", 1, i, j, LOG_PATH)
    lht_edata = list(hard = round(lht_hard$EstimationData, digits = 3), soft = round(lht_soft$EstimationData, digits = 3))
    lht_txt_path = create_txt(i,j,directory_path,"lht","B1")
    write.table(lht_edata$hard,lht_txt_path$hdata, row.names = FALSE)
    write.table(lht_edata$soft,lht_txt_path$sdata, row.names = FALSE)
  }
  
  # HBT_B2
  directory_path = "./output/DT_Bar_WSE/B2/"
  dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
  for(i in 2:log(GetGroupLength(length(ds)), base = 2)){
    ut_hard = run_wse_logged(ds, "B2", "ut", "h", 1, i, j, LOG_PATH)
    ut_soft = run_wse_logged(ds, "B2", "ut", "s", 1, i, j, LOG_PATH)
    ut_edata = list(hard = round(ut_hard$EstimationData, digits = 3), soft = round(ut_soft$EstimationData, digits = 3))
    ut_txt_path = create_txt(i,j,directory_path,"ut","B2")
    write.table(ut_edata$hard,ut_txt_path$hdata, row.names = FALSE)
    write.table(ut_edata$soft,ut_txt_path$sdata, row.names = FALSE)
    
    lut_hard = run_wse_logged(ds, "B2", "lut", "h", 1, i, j, LOG_PATH)
    lut_soft = run_wse_logged(ds, "B2", "lut", "s", 1, i, j, LOG_PATH)
    lut_edata = list(hard = round(lut_hard$EstimationData, digits = 3), soft = round(lut_soft$EstimationData, digits = 3))
    lut_txt_path = create_txt(i,j,directory_path,"lut","B2")
    write.table(lut_edata$hard,lut_txt_path$hdata, row.names = FALSE)
    write.table(lut_edata$soft,lut_txt_path$sdata, row.names = FALSE)
    
    lht_hard = run_wse_logged(ds, "B2", "lht", "h", 1, i, j, LOG_PATH)
    lht_soft = run_wse_logged(ds, "B2", "lht", "s", 1, i, j, LOG_PATH)
    lht_edata = list(hard = round(lht_hard$EstimationData, digits = 3), soft = round(lht_soft$EstimationData, digits = 3))
    lht_txt_path = create_txt(i,j,directory_path,"lht","B2")
    write.table(lht_edata$hard,lht_txt_path$hdata, row.names = FALSE)
    write.table(lht_edata$soft,lht_txt_path$sdata, row.names = FALSE)
  }
  
  # ===============================================================
  # HFreT (Freeman)
  # ===============================================================
  directory_path = "./output/DT_Fre_WSE/"
  dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
  for(i in 2:log(GetGroupLength(length(ds)), base = 2)){
    ut_hard = run_wse_logged(ds, "Fr", "ut", "h", 1, i, j, LOG_PATH)
    ut_soft = run_wse_logged(ds, "Fr", "ut", "s", 1, i, j, LOG_PATH)
    ut_edata = list(hard = round(ut_hard$EstimationData, digits = 3), soft = round(ut_soft$EstimationData, digits = 3))
    ut_txt_path = create_txt(i,j,directory_path,"ut","Fr")
    write.table(ut_edata$hard,ut_txt_path$hdata, row.names = FALSE)
    write.table(ut_edata$soft,ut_txt_path$sdata, row.names = FALSE)
    
    lut_hard = run_wse_logged(ds, "Fr", "lut", "h", 1, i, j, LOG_PATH)
    lut_soft = run_wse_logged(ds, "Fr", "lut", "s", 1, i, j, LOG_PATH)
    lut_edata = list(hard = round(lut_hard$EstimationData, digits = 3), soft = round(lut_soft$EstimationData, digits = 3))
    lut_txt_path = create_txt(i,j,directory_path,"lut","Fr")
    write.table(lut_edata$hard,lut_txt_path$hdata, row.names = FALSE)
    write.table(lut_edata$soft,lut_txt_path$sdata, row.names = FALSE)
    
    lht_hard = run_wse_logged(ds, "Fr", "lht", "h", 1, i, j, LOG_PATH)
    lht_soft = run_wse_logged(ds, "Fr", "lht", "s", 1, i, j, LOG_PATH)
    lht_edata = list(hard = round(lht_hard$EstimationData, digits = 3), soft = round(lht_soft$EstimationData, digits = 3))
    lht_txt_path = create_txt(i,j,directory_path,"lht","Fr")
    write.table(lht_edata$hard,lht_txt_path$hdata, row.names = FALSE)
    write.table(lht_edata$soft,lht_txt_path$sdata, row.names = FALSE)
  }
  
  # ===============================================================
  # HFitT (Fisz)
  # ===============================================================
  directory_path = "./output/DT_Fit_WSE/"
  dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
  for(i in 2:log(GetGroupLength(length(ds)), base = 2)){
    ut_hard = run_wse_logged(ds, "Fi", "ut", "h", 1, i, j, LOG_PATH)
    ut_soft = run_wse_logged(ds, "Fi", "ut", "s", 1, i, j, LOG_PATH)
    ut_edata = list(hard = round(ut_hard$EstimationData, digits = 3), soft = round(ut_soft$EstimationData, digits = 3))
    ut_txt_path = create_txt(i,j,directory_path,"ut","Fi")
    write.table(ut_edata$hard,ut_txt_path$hdata, row.names = FALSE)
    write.table(ut_edata$soft,ut_txt_path$sdata, row.names = FALSE)
    
    lut_hard = run_wse_logged(ds, "Fi", "lut", "h", 1, i, j, LOG_PATH)
    lut_soft = run_wse_logged(ds, "Fi", "lut", "s", 1, i, j, LOG_PATH)
    lut_edata = list(hard = round(lut_hard$EstimationData, digits = 3), soft = round(lut_soft$EstimationData, digits = 3))
    lut_txt_path = create_txt(i,j,directory_path,"lut","Fi")
    write.table(lut_edata$hard,lut_txt_path$hdata, row.names = FALSE)
    write.table(lut_edata$soft,lut_txt_path$sdata, row.names = FALSE)
    
    lht_hard = run_wse_logged(ds, "Fi", "lht", "h", 1, i, j, LOG_PATH)
    lht_soft = run_wse_logged(ds, "Fi", "lht", "s", 1, i, j, LOG_PATH)
    lht_edata = list(hard = round(lht_hard$EstimationData, digits = 3), soft = round(lht_soft$EstimationData, digits = 3))
    lht_txt_path = create_txt(i,j,directory_path,"lht","Fi")
    write.table(lht_edata$hard,lht_txt_path$hdata, row.names = FALSE)
    write.table(lht_edata$soft,lht_txt_path$sdata, row.names = FALSE)
  }
}

