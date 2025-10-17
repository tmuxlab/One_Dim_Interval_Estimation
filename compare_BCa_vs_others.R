# ============================================================
# BCa vs PercAns / tQuick 比較（幅・位置・重なり・被覆率）
# 依存：あなたの規約で出力された
#   - Main1:   ./output/<GroupDir>/<FileName?>/DS*-k<dt>-<thrMode>-<thrName>(J=*).txt
#   - PE:      ./PE-output/<GroupDir>/<FileName?>/DS*-PE-<BSMethod>-k<dt>-...dat
#   - IE(Perc/tQuick): ./IE-output/<GroupDir>/<FileName?>/DS*-IEEachPnt-<PercAns|tQuick>-...csv
#   - JK:      ./output/<GroupDir>/<FileName?>/JK-intensity/DS*-k<dt>-...-JK-<i>.txt（既存 or 生成済み）
# ============================================================

suppressPackageStartupMessages({
  library(stats); library(utils)
})

# --- ユーティリティ（読み込み） ---
.safe_path <- function(...) { p <- c(...); p <- p[nzchar(p)]; do.call(file.path, as.list(p)) }

read_numeric_col <- function(path){
  L <- readLines(path, warn=FALSE)
  L <- L[!grepl("^\\s*$", L)]; L <- L[!grepl("^\\s*#", L)]; L <- trimws(L)
  tok <- sub("^([^\\t ,;]+).*","\\1", L)
  v <- suppressWarnings(as.numeric(tok)); v <- v[is.finite(v)]
  if (!length(v)) stop(sprintf("数値を読めませんでした: %s", path))
  v
}

read_numeric_matrix <- function(path){
  raw <- suppressWarnings(read.table(path, header=FALSE, sep="", quote="", comment.char="",
                                     check.names=FALSE, stringsAsFactors=FALSE, fill=TRUE))
  for (j in seq_along(raw)) raw[[j]] <- suppressWarnings(as.numeric(raw[[j]]))
  M <- as.matrix(raw)
  keep_col <- colSums(is.finite(M))>0; keep_row <- rowSums(is.finite(M))>0
  M <- M[keep_row, keep_col, drop=FALSE]
  if (!is.matrix(M) || any(dim(M)==0)) stop(sprintf("数値行列を構成できませんでした: %s", path))
  M
}

read_jk_matrix <- function(dir_path, DS, dt_code, thrMode, thrName, J){
  # 例: DS1-kB2-h-ut(J=2)-JK-<i>.txt
  pat <- sprintf("^DS%s-k%s-%s-%s\\(J=%s\\)-JK-(\\d+)\\.txt$", DS, dt_code, thrMode, thrName, J)
  files <- list.files(dir_path, pattern = pat, full.names = TRUE)
  if (!length(files)) return(NULL)
  rx <- regexec(pat, basename(files)); m <- regmatches(basename(files), rx)
  idx <- sapply(m, function(mm) as.integer(mm[2])); o <- order(idx); files <- files[o]
  v0 <- as.numeric(read.table(files[1])[[1]]); Tlen <- length(v0)
  JK <- matrix(NA_real_, nrow = length(files), ncol = Tlen)
  for (i in seq_along(files)) JK[i,] <- as.numeric(read.table(files[i])[[1]])
  JK
}

# --- BCa 分位（全時点） ---
bca_quantile_vector <- function(boot_vec, theta_hat, alpha=0.05, a=0){
  v <- boot_vec[is.finite(boot_vec)]
  if (length(v) < 2 || !is.finite(theta_hat)) return(c(NA_real_, NA_real_))
  prop_lt <- mean(v < theta_hat)
  prop_eq <- mean(v == theta_hat)
  z0 <- qnorm(prop_lt + 0.5*prop_eq)
  zL <- qnorm(alpha/2); zU <- qnorm(1 - alpha/2)
  adj <- function(z) pnorm( z0 + (z0+z) / (1 - a*(z0+z)) )
  pL <- min(max(adj(zL), 0), 1); pU <- min(max(adj(zU), 0), 1)
  as.numeric(quantile(v, probs=c(pL,pU), type=7, names=FALSE, na.rm=TRUE))
}

bca_ci_from_PE <- function(PE_path, lambda_hat, JK=NULL, alpha=0.05){
  M <- read_numeric_matrix(PE_path)   # BxT or TxB を許容
  Tlen <- length(lambda_hat)
  if (ncol(M)==Tlen) { B <- nrow(M) } else if (nrow(M)==Tlen) { M <- t(M); B <- nrow(M) } else {
    stop(sprintf("PE形状が T=%d と不一致: %dx%d", Tlen, nrow(M), ncol(M)))
  }
  # a（加速度）
  if (!is.null(JK)) {
    accel_from_jackknife_vec <- function(v){ th<-mean(v); u<-th - v; den<-6*(sum(u^2))^(3/2); if(!is.finite(den)||den==0) 0 else sum(u^3)/den }
    a <- apply(JK, 2, accel_from_jackknife_vec)  # 長さT
  } else {
    a <- rep(0, Tlen)
  }
  lo <- hi <- rep(NA_real_, Tlen)
  for (t in seq_len(Tlen)){
    q <- bca_quantile_vector(M[,t], lambda_hat[t], alpha=alpha, a=a[t])
    lo[t] <- q[1]; hi[t] <- q[2]
  }
  data.frame(lo=lo, hi=hi)
}

# --- 比較（幅・位置・重なり・被覆率） ---
ci_compare_two <- function(ci_ref, ci_new, lambda_hat, alpha=0.05){
  stopifnot(all(c("lo","hi") %in% names(ci_ref)), all(c("lo","hi") %in% names(ci_new)))
  width_ref <- pmax(0, ci_ref$hi - ci_ref$lo)
  width_new <- pmax(0, ci_new$hi - ci_new$lo)
  mid_ref   <- (ci_ref$hi + ci_ref$lo)/2
  mid_new   <- (ci_new$hi + ci_new$lo)/2
  
  overlap <- pmax(0, pmin(ci_ref$hi, ci_new$hi) - pmax(ci_ref$lo, ci_new$lo))
  union   <- pmax(width_ref + width_new - overlap, 1e-12)
  jacc    <- overlap / union
  
  inside  <- (lambda_hat >= ci_new$lo) & (lambda_hat <= ci_new$hi)
  bt      <- binom.test(sum(inside, na.rm=TRUE), length(lambda_hat), p=1-alpha)
  
  list(
    mean_width_ref = mean(width_ref, na.rm=TRUE),
    mean_width_new = mean(width_new, na.rm=TRUE),
    mean_width_diff = mean(width_new - width_ref, na.rm=TRUE),
    prop_new_wider  = mean(width_new > width_ref, na.rm=TRUE),
    mean_mid_absdiff= mean(abs(mid_new - mid_ref), na.rm=TRUE),
    rmse_mid_diff   = sqrt(mean((mid_new - mid_ref)^2, na.rm=TRUE)),
    mean_overlap_ratio = mean(jacc, na.rm=TRUE),
    coverage_inside = sum(inside, na.rm=TRUE),
    coverage_prop   = mean(inside, na.rm=TRUE),
    coverage_p      = bt$p.value
  )
}

# --- 主要関数：BCa vs PercAns / tQuick ---
compare_BCa_vs_others <- function(
    DS, GroupDir, FileName = "",
    BSMethod, dt_code, thrMode, thrName, J,
    cos_tag = "COS",
    alpha = 0.05,
    jk_dir = NULL,                 # 例: "./output/DT_Bar_WSE/B2/JK-intensity"
    make_plot = FALSE              # TRUEで重ね描きPNG出力
){
  # compare_BCa_vs_others(...) の冒頭で
  if (missing(BSMethod) || is.null(BSMethod)) BSMethod <- cos_tag
  if (!identical(BSMethod, cos_tag)) {
    stop("BSMethod (PEのタグ) と cos_tag (IEのタグ) が不一致です：",
         "BCaの基礎PEと比較先IEが別COSになってしまいます。BSMethod=cos_tagにしてください。")
  }
  
  # パス
  base_out   <- .safe_path("./output",   GroupDir, FileName)
  base_pe    <- .safe_path("./PE-output",GroupDir, FileName)
  base_ie    <- .safe_path("./IE-output",GroupDir, FileName)
  
  kcode_eff  <- paste0("k", dt_code)
  main1_path <- .safe_path(base_out, sprintf("DS%s-%s-%s-%s(J=%s).txt", DS, kcode_eff, thrMode, thrName, J))
  pe_I_path  <- .safe_path(base_pe, sprintf("DS%s-PE-%s-k%s-%s-%s(J=%s).dat", DS, BSMethod, dt_code, thrMode, thrName, J))
  # 既存CI
  f_Perc <- .safe_path(base_ie, sprintf("DS%s-IEEachPnt-PercAns-%s-%s-%s-%s(J=%s).csv", DS, cos_tag, kcode_eff, thrMode, thrName, J))
  f_tQ   <- .safe_path(base_ie, sprintf("DS%s-IEEachPnt-tQuick-%s-%s-%s-%s(J=%s).csv",  DS, cos_tag, kcode_eff, thrMode, thrName, J))
  if (!file.exists(f_Perc)) stop("PercAns CI not found: ", f_Perc)
  if (!file.exists(f_tQ))   stop("tQuick CI not found: ",   f_tQ)
  
  lambda_hat <- read_numeric_col(main1_path)
  ci_perc <- read.csv(f_Perc, check.names=FALSE)
  ci_tq   <- read.csv(f_tQ,   check.names=FALSE)
  
  # JK 読み込み（なければ加速度 a=0 で実行）
  if (is.null(jk_dir)) jk_dir <- .safe_path(base_out, "JK-intensity")
  JK <- if (dir.exists(jk_dir)) read_jk_matrix(jk_dir, DS, dt_code, thrMode, thrName, J) else NULL
  
  # BCa CI（全時点）を PE から生成
  ci_bca <- bca_ci_from_PE(pe_I_path, lambda_hat, JK, alpha=alpha)
  
  # 比較：BCa vs PercAns / tQuick
  cmp_vs_perc <- ci_compare_two(ci_ref = ci_perc, ci_new = ci_bca, lambda_hat = lambda_hat, alpha = alpha)
  cmp_vs_tq   <- ci_compare_two(ci_ref = ci_tq,   ci_new = ci_bca, lambda_hat = lambda_hat, alpha = alpha)
  
  # まとめテーブル（1行）
  out_row <- data.frame(
    DS=DS, GroupDir=GroupDir, FileName=FileName, COS=cos_tag, J=J,
    BSMethod=BSMethod, dt=dt_code, thrMode=thrMode, thrName=thrName,
    # PercAns基準
    BCA_vs_Perc_mean_width_ref = cmp_vs_perc$mean_width_ref,
    BCA_vs_Perc_mean_width_new = cmp_vs_perc$mean_width_new,
    BCA_vs_Perc_mean_width_diff= cmp_vs_perc$mean_width_diff,
    BCA_vs_Perc_prop_new_wider = cmp_vs_perc$prop_new_wider,
    BCA_vs_Perc_rmse_mid_diff  = cmp_vs_perc$rmse_mid_diff,
    BCA_vs_Perc_overlap_ratio  = cmp_vs_perc$mean_overlap_ratio,
    BCA_coverage_inside        = cmp_vs_perc$coverage_inside,
    BCA_coverage_prop          = cmp_vs_perc$coverage_prop,
    BCA_coverage_p             = cmp_vs_perc$coverage_p,
    # tQuick基準
    BCA_vs_tQ_mean_width_ref   = cmp_vs_tq$mean_width_ref,
    BCA_vs_tQ_mean_width_new   = cmp_vs_tq$mean_width_new,
    BCA_vs_tQ_mean_width_diff  = cmp_vs_tq$mean_width_diff,
    BCA_vs_tQ_prop_new_wider   = cmp_vs_tq$prop_new_wider,
    BCA_vs_tQ_rmse_mid_diff    = cmp_vs_tq$rmse_mid_diff,
    BCA_vs_tQ_overlap_ratio    = cmp_vs_tq$mean_overlap_ratio,
    stringsAsFactors = FALSE
  )
  
  # 保存
  out_dir <- base_ie
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  out_csv <- .safe_path(out_dir, sprintf("DS%s-BCACOMP-%s-%s-%s-%s(J=%s).csv", DS, cos_tag, kcode_eff, thrMode, thrName, J))
  write.csv(out_row, out_csv, row.names=FALSE)
  message("Saved: ", out_csv)
  
  # 可視化（任意）：λ̂ と3種CI(BCA/Perc/tQuick)の重ね描き
  if (isTRUE(make_plot)){
    pngfile <- .safe_path(out_dir, sprintf("DS%s-BCACOMP-PLOT-%s-%s-%s-%s(J=%s).png", DS, cos_tag, kcode_eff, thrMode, thrName, J))
    png(pngfile, width=1400, height=900, res=150)
    Tlen <- length(lambda_hat); plot(lambda_hat, type="l", col=2, lwd=2, xlab="t", ylab="lambda",
                                     main=sprintf("BCa vs PercAns / tQuick (%s, J=%s)", cos_tag, J))
    lines(ci_perc$lo, lty=2); lines(ci_perc$hi, lty=2)
    lines(ci_tq$lo,   lty=3, col=4); lines(ci_tq$hi, lty=3, col=4)
    lines(ci_bca$lo,  lty=1, col=1); lines(ci_bca$hi, lty=1, col=1)
    legend("topright",
           c("lambda_hat","PercAns","tQuick","BCa"),
           lty=c(1,2,3,1), col=c(2,1,4,1), lwd=c(2,1,1,1), bty="n")
    dev.off(); message("Saved plot: ", pngfile)
  }
  
  invisible(list(csv=out_csv, row=out_row, ci_bca=ci_bca))
}
