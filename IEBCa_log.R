## =========================================
##  BCa一括実行：JKをオンザフライ作成＋Jループ対応
##  前提：あなたの環境で Wse() / create_txt() / GetGroupLength() が使える
##        run_wse_logged() もあなたの定義を再利用（別ログCSVへ）
## =========================================

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

# Load Hal wavelet estimation module
WSE_Path = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/src/WaveletShrinkageEstimation.R")
print("Load Hal wavelet estimation module")
source(WSE_Path)

WT_Path = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/src/WaveletTransform.R")
print("Load Hal wavelet Transformation module")
source(WT_Path)

suppressPackageStartupMessages(library(tidyverse))

## -------- パス/命名ユーティリティ --------
.safe_path <- function(...) { p <- c(...); p <- p[nzchar(p)]; do.call(file.path, as.list(p)) }

# Main1 出力ファイル（観測データ全量での強度推定）
main1_out_path <- function(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, J) {
        .safe_path("./output", GroupDir, FileName,
                   sprintf("DS%s-k%s-%s-%s(J=%s).txt", DS, dt_code, thresholdMode, thresholdName, J))
}

# PEファイル（あなたの実ファイル名規約に合わせる：強度=k{dt} / 累積=r{dt}）
pe_paths <- function(DS, GroupDir, FileName, BSMethod, dt_code, thresholdMode, thresholdName, J) {
        ped <- .safe_path("./PE-output", GroupDir, FileName)
        list(
                intensity = .safe_path(ped, sprintf("DS%s-PE-%s-k%s-%s-%s(J=%s).dat", DS, BSMethod, dt_code, thresholdMode, thresholdName, J)),
                cumulative = .safe_path(ped, sprintf("DS%s-PE-%s-r%s-%s-%s(J=%s).dat", DS, BSMethod, dt_code, thresholdMode, thresholdName, J))
        )
}

# IE出力タグ
tag_str <- function(BSMethod, dt_code, thresholdMode, thresholdName, J) {
        sprintf("%s-k%s-%s-%s(J=%s)", BSMethod, dt_code, thresholdMode, thresholdName, J)
}

## -------- あなたの run_wse_logged() を活用して JK を作る --------
# dt_code → Wse() の transform へ（H は "none"、A1/A2/A3/B1/B2 はそのまま）
dt_to_transform <- function(dt_code) {
        if (dt_code %in% c("H", "NDT", "NONE", "none", "h")) return("none")
        dt_code
}

# JKログは本線と別ファイルに
LOG_PATH_JK <- "./output/metrics/wse_jk_metrics.csv"

# 1回の WSE 推定（ds: 観測カウントベクトル）
fit_wse_once <- function(ds, dt_code, thresholdName, thresholdMode, J, dataset_id) {
        transform <- dt_to_transform(dt_code)
        res <- run_wse_logged(ds, transform, thresholdName, thresholdMode, 1, J, dataset_id, LOG_PATH_JK)
        as.numeric(res$EstimationData)
}

# JK をオンザフライで生成（保存もする）
build_jk_intensity <- function(
                counts_vec, DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, J,
                save_dir = NULL, save_files = TRUE
){
        n <- length(counts_vec)
        JKdir <- if (is.null(save_dir)) .safe_path("./output", GroupDir, FileName, "JK-intensity") else save_dir
        if (save_files) dir.create(JKdir, recursive = TRUE, showWarnings = FALSE)
        
        JK <- NULL
        for (i in seq_len(n)) {
                ds_loo <- counts_vec[-i]
                vhat <- fit_wse_once(ds_loo, dt_code, thresholdName, thresholdMode, J, dataset_id = DS)
                if (i == 1) JK <- matrix(NA_real_, nrow = n, ncol = length(vhat))
                JK[i, ] <- vhat
                
                if (save_files) {
                        fn <- sprintf("DS%s-k%s-%s-%s(J=%s)-JK-%s.txt", DS, dt_code, thresholdMode, thresholdName, J, i)
                        write.table(vhat, file = .safe_path(JKdir, fn), row.names = FALSE, col.names = FALSE, quote = FALSE)
                }
        }
        JK
}

# 強度JK → 累積JK（行=LOO、列=時点）
jk_cum_from_intensity <- function(JK) t(apply(JK, 1, cumsum))

## ===== (1) 計測ラッパ：BS_intervals と同等の拡張版 =====
.get_rss_mb <- function() {
        sys <- tolower(Sys.info()[["sysname"]]); pid <- Sys.getpid()
        parse_mb <- function(xkb){ v <- suppressWarnings(as.numeric(xkb)); if(!is.finite(v)) return(NA_real_); v/1024 }
        if (grepl("linux|darwin", sys)) {
                out <- suppressWarnings(tryCatch(system(sprintf("ps -o rss= -p %d", pid), intern=TRUE), error=function(e) NA_character_))
                if (length(out) && !is.na(out[1])) return(parse_mb(out[1]))
                if (grepl("linux", sys)) {
                        p <- sprintf("/proc/%d/status", pid); if (file.exists(p)) {
                                ln <- tryCatch(readLines(p, warn=FALSE), error=function(e) character(0))
                                vm <- ln[grepl("^VmRSS:", ln)]; if (length(vm)) return(parse_mb(sub(".*?([0-9]+).*","\\1", vm[1])))
                        }
                }
        } else if (grepl("windows", sys)) {
                cmd <- sprintf('powershell -NoProfile -Command "(Get-Process -Id %d).WorkingSet64"', pid)
                out <- suppressWarnings(tryCatch(system(cmd, intern=TRUE), error=function(e) NA_character_))
                if (length(out) && !is.na(out[1])) { b <- suppressWarnings(as.numeric(out[1])); if (is.finite(b)) return(b/(1024^2)) }
        }
        NA_real_
}
.mem_used_mb <- function() {
        if (requireNamespace("pryr", quietly=TRUE)) {
                v <- as.numeric(pryr::mem_used())/(1024^2); if (!is.finite(v)) return(NA_real_); return(v)
        }
        if (.Platform$OS.type=="windows") { v <- suppressWarnings(tryCatch(utils::memory.size(), error=function(e) NA_real_)); if (!is.finite(v)) v <- NA_real_; return(v) }
        NA_real_
}
with_perf_log_extended <- function(tag, meta=list(), code, log_dir="./logs", log_name="perf_log_extended.csv"){
        code_sub <- substitute(code)
        rss0 <- .get_rss_mb(); mu0 <- .mem_used_mb(); t0 <- proc.time()
        result <- eval(code_sub, envir = parent.frame())
        pt <- proc.time() - t0; rss1 <- .get_rss_mb(); mu1 <- .mem_used_mb()
        nf <- function(x){ x <- suppressWarnings(as.numeric(x)); ifelse(is.finite(x), x, NA_real_) }
        row <- data.frame(
                dataset = if (!is.null(meta$dataset)) meta$dataset else NA,
                J = if (!is.null(meta$J)) meta$J else NA,
                FileName = if (!is.null(meta$FileName)) meta$FileName else "",
                BSMethod = if (!is.null(meta$BSMethod)) meta$BSMethod else "",
                dt = if (!is.null(meta$dt)) meta$dt else "",
                thresholdMode = if (!is.null(meta$thresholdMode)) meta$thresholdMode else "",
                thresholdName = if (!is.null(meta$thresholdName)) meta$thresholdName else "",
                timestamp = as.character(Sys.time()), tag = as.character(tag),
                user_sec = nf(pt[["user.self"]]), system_sec = nf(pt[["sys.self"]]), elapsed_sec = nf(pt[["elapsed"]]),
                rss_mb_before = nf(rss0), rss_mb_after = nf(rss1),
                rss_mb_delta = if (is.finite(rss0) && is.finite(rss1)) rss1-rss0 else NA_real_,
                mem_used_mb_before = nf(mu0), mem_used_mb_after = nf(mu1),
                mem_used_mb_delta = if (is.finite(mu0) && is.finite(mu1)) mu1-mu0 else NA_real_,
                stringsAsFactors = FALSE
        )
        dir.create(log_dir, recursive=TRUE, showWarnings=FALSE)
        log_path <- file.path(log_dir, log_name)
        utils::write.table(row, file=log_path, sep=",", quote=TRUE, row.names=FALSE,
                           col.names=!file.exists(log_path), append=TRUE)
        invisible(result)
}

## ===== (2) 読み込みヘルパ（堅牢） =====
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

## ===== (3) BCa 本体＋計測・幅ファイルの追加 =====
IE_BCa_from_PE2 <- function(
                DS, GroupDir, FileName,
                BSMethod, dt_code, thresholdMode, thresholdName, J,
                obs_estimate_path,
                jackknife_intensity = NULL, jackknife_cum = NULL,
                alpha = 0.05, tiny = 1e-8
){
        # パス
        pe_dir <- file.path("./PE-output", GroupDir, if (nzchar(FileName)) FileName else "")
        peI <- file.path(pe_dir, sprintf("DS%s-PE-%s-k%s-%s-%s(J=%s).dat", DS, BSMethod, dt_code, thresholdMode, thresholdName, J))
        peC <- file.path(pe_dir, sprintf("DS%s-PE-%s-r%s-%s-%s(J=%s).dat", DS, BSMethod, dt_code, thresholdMode, thresholdName, J))
        out_dir <- file.path("./IE-output", GroupDir, if (nzchar(FileName)) FileName else "")
        dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
        tag <- sprintf("%s-k%s-%s-%s(J=%s)", BSMethod, dt_code, thresholdMode, thresholdName, J)
        
        meta <- list(dataset=DS, J=J, FileName=FileName, BSMethod=BSMethod, dt=dt_code,
                     thresholdMode=thresholdMode, thresholdName=thresholdName)
        
        with_perf_log_extended("IE_BCa", meta, {
                # 読み込み
                boots_I   <- read_numeric_matrix(peI)
                boots_Cum <- read_numeric_matrix(peC)
                est_vec   <- read_numeric_col(obs_estimate_path)
                T_len <- length(est_vec)
                to_BxT <- function(M){ if(ncol(M)==T_len) M else if(nrow(M)==T_len) t(M) else stop(sprintf("PE 形状が T=%d と不一致: %dx%d", T_len, nrow(M), ncol(M))) }
                boots_I   <- to_BxT(boots_I); boots_Cum <- to_BxT(boots_Cum)
                
                # a（加速度）
                accel_from_jackknife_vec <- function(v){ th<-mean(v); u<-th - v; den<-6*(sum(u^2))^(3/2); if(!is.finite(den)||den==0) 0 else sum(u^3)/den }
                accel_from_jackknife_mat <- function(JK) apply(JK, 2, accel_from_jackknife_vec)
                a_I   <- if (is.null(jackknife_intensity)) 0 else accel_from_jackknife_mat(jackknife_intensity)
                a_Cum <- if (is.null(jackknife_cum))       0 else accel_from_jackknife_mat(jackknife_cum)
                
                # BCa 分位
                BCa_quant <- function(boot_vec, theta_hat, alpha, a=0){
                        v <- boot_vec[is.finite(boot_vec)]; if (length(v) < 2) return(c(NA_real_, NA_real_))
                        prop_lt <- mean(v < theta_hat); prop_eq <- mean(v == theta_hat)
                        z0 <- qnorm(prop_lt + 0.5*prop_eq); zL <- qnorm(alpha/2); zU <- qnorm(1 - alpha/2)
                        adj <- function(z) pnorm( z0 + (z0+z) / (1 - a*(z0+z)) )
                        pL <- min(max(adj(zL), 0), 1); pU <- min(max(adj(zU), 0), 1)
                        as.numeric(quantile(v, probs=c(pL,pU), type=7, names=FALSE, na.rm=TRUE))
                }
                
                # 並べ替え保存（IEEach）
                ord_I   <- t(apply(boots_I,   2, sort))
                ord_Cum <- t(apply(boots_Cum, 2, sort))
                f_IEEach_I   <- file.path(out_dir, sprintf("DS%s-IEEachPnt-%s.csv",  DS, tag))
                f_IEEach_Cum <- file.path(out_dir, sprintf("DS%s-IEEachPntC-%s.csv", DS, tag))
                write.csv(as.data.frame(ord_I),   f_IEEach_I,   quote=FALSE, row.names=FALSE)
                write.csv(as.data.frame(ord_Cum), f_IEEach_Cum, quote=FALSE, row.names=FALSE)
                
                # 4点（25/50/75/100%時点）要約 + BCa CI
                if (!requireNamespace("e1071", quietly=TRUE)) stop("e1071 を install.packages('e1071') で導入してください。")
                idx <- pmax(1, pmin(T_len, floor(c(0.25, 0.5, 0.75, 1.0) * T_len)))
                MomIE  <- matrix(NA_real_, 4, 7); MomIEc <- matrix(NA_real_, 4, 7)
                # 幅（UL-LL）：四分位点／全時点
                w_ie_q  <- numeric(length(idx)); w_iec_q <- numeric(length(idx))
                w_ie_all  <- numeric(T_len);     w_iec_all <- numeric(T_len)
                
                for (k in seq_along(idx)) {
                        j <- idx[k]
                        vI   <- boots_I[,j];   vI   <- vI  [vI   != tiny & is.finite(vI)]
                        vCum <- boots_Cum[,j]; vCum <- vCum[vCum != tiny & is.finite(vCum)]
                        thI   <- est_vec[j];   thCum <- sum(est_vec[seq_len(j)])
                        aI   <- if (length(a_I)>1)   a_I[j]   else a_I
                        aCum <- if (length(a_Cum)>1) a_Cum[j] else a_Cum
                        ciI   <- BCa_quant(vI,   thI,   alpha, aI)
                        ciCum <- BCa_quant(vCum, thCum, alpha, aCum)
                        
                        MomIE [k,] <- c(round(mean(vI),4),   round(median(vI),4),   round(var(vI),4),
                                        round(e1071::skewness(vI),4), round(e1071::kurtosis(vI)+3,4),
                                        round(ciI[1],4),   round(ciI[2],4))
                        MomIEc[k,] <- c(round(mean(vCum),4), round(median(vCum),4), round(var(vCum),4),
                                        round(e1071::skewness(vCum),4), round(e1071::kurtosis(vCum)+3,4),
                                        round(ciCum[1],4), round(ciCum[2],4))
                        
                        w_ie_q[k]  <- as.numeric(ciI[2]  - ciI[1])
                        w_iec_q[k] <- as.numeric(ciCum[2] - ciCum[1])
                }
                
                # 全時点の幅（UL-LL）
                probs <- c(0.025, 0.975)
                for (j in seq_len(T_len)) {
                        vv  <- boots_I[,j];   vv  <- vv [vv  != tiny & is.finite(vv)]
                        vvc <- boots_Cum[,j]; vvc <- vvc[vvc != tiny & is.finite(vvc)]
                        qv  <- if (length(vv))  quantile(vv,  probs=probs, type=7, names=FALSE) else c(NA_real_, NA_real_)
                        qvc <- if (length(vvc)) quantile(vvc, probs=probs, type=7, names=FALSE) else c(NA_real_, NA_real_)
                        w_ie_all[j]  <- as.numeric(qv[2]  - qv[1])
                        w_iec_all[j] <- as.numeric(qvc[2] - qvc[1])
                }
                
                # IE（MomIE/MomIEc）
                f_MomIE_BCa  <- file.path(out_dir, sprintf("DS%s-IEBCa-%s.dat",  DS, tag))
                f_MomIEc_BCa <- file.path(out_dir, sprintf("DS%s-IEcBCa-%s.dat", DS, tag))
                write.table(MomIE,  f_MomIE_BCa,  quote=FALSE, col.names=FALSE, row.names=FALSE)
                write.table(MomIEc, f_MomIEc_BCa, quote=FALSE, col.names=FALSE, row.names=FALSE)
                
                # 幅の表・要約（IE_intervals と同趣旨のファイル群）
                pcnts <- c(0.25, 0.5, 0.75, 1.0)
                ie_width_tbl  <- cbind(TimeFrac = pcnts, CI_Width = round(w_ie_q,  6))
                iec_width_tbl <- cbind(TimeFrac = pcnts, CI_Width = round(w_iec_q, 6))
                
                f_IEw        <- file.path(out_dir, sprintf("DS%s-IEBCawidth-%s.dat",         DS, tag))
                f_IEcw       <- file.path(out_dir, sprintf("DS%s-IEcBCawidth-%s.dat",        DS, tag))
                f_IEw_sum    <- file.path(out_dir, sprintf("DS%s-IEBCawidth-summary-%s.dat", DS, tag))
                f_IEcw_sum   <- file.path(out_dir, sprintf("DS%s-IEcBCawidth-summary-%s.dat",DS, tag))
                f_IEw_allsum <- file.path(out_dir, sprintf("DS%s-IEBCawidth-all-summary-%s.dat",  DS, tag))
                f_IEcw_allsum<- file.path(out_dir, sprintf("DS%s-IEcBCawidth-all-summary-%s.dat", DS, tag))
                
                write.table(ie_width_tbl,  f_IEw,  quote=FALSE, col.names=TRUE,  row.names=FALSE)
                write.table(iec_width_tbl, f_IEcw, quote=FALSE, col.names=TRUE,  row.names=FALSE)
                
                ie_w_mean <- round(mean(w_ie_q,  na.rm=TRUE), 6); ie_w_var <- round(var(w_ie_q,   na.rm=TRUE), 6); ie_w_sd <- round(sd(w_ie_q,   na.rm=TRUE), 6)
                ic_w_mean <- round(mean(w_iec_q, na.rm=TRUE), 6); ic_w_var <- round(var(w_iec_q,  na.rm=TRUE), 6); ic_w_sd <- round(sd(w_iec_q,  na.rm=TRUE), 6)
                write.table(matrix(c(ie_w_mean,  ie_w_var,  ie_w_sd),  nrow=1), f_IEw_sum,  quote=FALSE, col.names=FALSE, row.names=FALSE)
                write.table(matrix(c(ic_w_mean,  ic_w_var,  ic_w_sd),  nrow=1), f_IEcw_sum, quote=FALSE, col.names=FALSE, row.names=FALSE)
                
                # 全時点の幅の要約
                ie_all_mean <- round(mean(w_ie_all,  na.rm=TRUE), 6); ie_all_var <- round(var(w_ie_all,   na.rm=TRUE), 6); ie_all_sd <- round(sd(w_ie_all,   na.rm=TRUE), 6)
                ic_all_mean <- round(mean(w_iec_all, na.rm=TRUE), 6); ic_all_var <- round(var(w_iec_all,  na.rm=TRUE), 6); ic_all_sd <- round(sd(w_iec_all,  na.rm=TRUE), 6)
                write.table(matrix(c(ie_all_mean,  ie_all_var,  ie_all_sd),  nrow=1), f_IEw_allsum,  quote=FALSE, col.names=FALSE, row.names=FALSE)
                write.table(matrix(c(ic_all_mean,  ic_all_var,  ic_all_sd),  nrow=1), f_IEcw_allsum, quote=FALSE, col.names=FALSE, row.names=FALSE)
        })
        
        invisible(list(out_dir=out_dir))
}

## -------- バッチ：J を回しつつ、JKがなければ作成してBCa実行 --------
run_BCa_batch <- function(
                DS, GroupDir, FileName = "",
                BSMethod, dt_code, thresholdMode, thresholdName,
                J_list,
                counts_path,                        # 観測データ（JK生成に必須）
                jk_file_dir = NULL, jk_regex_tpl = NULL,   # 既存JKを読むとき（例: "^DS1-kH-h-ldt\\(J={J}\\)-JK-(\\d+)\\.txt$")
                alpha = 0.05
){
        # 観測カウント（2列目が個数）
        raw <- read.table(counts_path, header = FALSE, check.names = FALSE)
        counts <- if (ncol(raw) >= 2) as.numeric(raw[[2]]) else as.numeric(raw[[1]])
        
        results <- list()
        for (J in J_list) {
                message(sprintf("=== J=%s ===", J))
                
                # Main1 推定ファイル（全データ）
                obs_est <- main1_out_path(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, J)
                
                # JK を用意（存在すれば読み込み・無ければオンザフライ作成＆保存）
                JK <- NULL
                if (!is.null(jk_file_dir) && !is.null(jk_regex_tpl) && dir.exists(jk_file_dir)) {
                        # 変更前（私の以前の案）
                        # regex_J <- gsub("\\{J\\}", J, jk_regex_tpl, fixed = TRUE)
                        
                        # 変更後（こちらにしてください）
                        regex_J <- gsub("{J}", as.character(J), jk_regex_tpl, fixed = TRUE)
                        
                        files <- list.files(jk_file_dir, pattern = regex_J, full.names = TRUE)
                        if (length(files)) {
                                # 読み込み（順序を正に）
                                rx <- regexec(regex_J, basename(files)); m <- regmatches(basename(files), rx)
                                idx <- sapply(m, function(mm) as.integer(mm[2])); o <- order(idx); files <- files[o]
                                # 最初で T 確認
                                v0 <- as.numeric(read.table(files[1])[[1]]); Tlen <- length(v0)
                                JK <- matrix(NA_real_, nrow = length(files), ncol = Tlen)
                                for (i in seq_along(files)) JK[i,] <- as.numeric(read.table(files[i])[[1]])
                        }
                }
                if (is.null(JK)) {
                        JK <- build_jk_intensity(
                                counts_vec = counts,
                                DS = DS, GroupDir = GroupDir, FileName = FileName,
                                dt_code = dt_code, thresholdMode = thresholdMode, thresholdName = thresholdName, J = J,
                                save_dir = .safe_path("./output", GroupDir, FileName, "JK-intensity"),
                                save_files = TRUE
                        )
                }
                JKc <- jk_cum_from_intensity(JK)
                
                # IE（BCa）
                res <- IE_BCa_from_PE2(
                        DS = DS, GroupDir = GroupDir, FileName = FileName,
                        BSMethod = BSMethod, dt_code = dt_code, thresholdMode = thresholdMode, thresholdName = thresholdName, J = J,
                        obs_estimate_path = obs_est,
                        jackknife_intensity = JK, jackknife_cum = JKc,
                        alpha = alpha
                )
                results[[as.character(J)]] <- list(J=J, IE=res, JK=JK)
        }
        invisible(results)
}

out <- run_BCa_batch(
  DS = 1,
  GroupDir = "DT_Bar_WSE",
  FileName = "B2",
  BSMethod = "ThiInt",
  dt_code = "B2",
  thresholdMode = "h",
  thresholdName = "ut",
  J_list = 2:5,
  counts_path = "./DS/DS1.txt",   # 観測カウント（2列目）
  jk_file_dir = "./output/DT_Bar_WSE/JK-intensity",
  jk_regex_tpl = "^DS1-kB2-h-ut\\(J={J}\\)-JK-(\\d+)\\.txt$",  # 既存があれば読込
  alpha = 0.05
)
