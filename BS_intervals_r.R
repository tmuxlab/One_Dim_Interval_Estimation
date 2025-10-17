# e1071 を明示して使う（他パッケージにマスクされないように :: を併用）
library(e1071)
# moments を既に読み込んでいる場合、不要なら外す（任意）
# if ("package:moments" %in% search()) detach("package:moments", unload=TRUE)
# ==== 計測ヘルパ =============================================

# ==== RSS (MB) を可能な限り堅牢に取得（非有限は NA に） ====
.get_rss_mb <- function() {
        sys <- tolower(Sys.info()[["sysname"]])
        pid <- Sys.getpid()
        
        parse_numeric_mb <- function(x_kb) {
                val <- suppressWarnings(as.numeric(x_kb))
                if (!is.finite(val)) return(NA_real_)
                val / 1024.0
        }
        
        if (grepl("linux", sys)) {
                # /proc/<pid>/status から VmRSS を読む
                path <- sprintf("/proc/%d/status", pid)
                if (file.exists(path)) {
                        ln <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
                        vmrss <- ln[grepl("^VmRSS:", ln)]
                        if (length(vmrss)) {
                                kb <- sub(".*?([0-9]+)\\s*kB.*", "\\1", vmrss[1])
                                return(parse_numeric_mb(kb))
                        }
                }
                # フォールバック: ps
                out <- suppressWarnings(tryCatch(system(sprintf("ps -o rss= -p %d", pid), intern = TRUE), error=function(e) NA_character_))
                if (length(out) && !is.na(out[1])) return(parse_numeric_mb(out[1]))
                return(NA_real_)
        }
        
        if (grepl("darwin", sys)) { # macOS
                out <- suppressWarnings(tryCatch(system(sprintf("ps -o rss= -p %d", pid), intern = TRUE), error=function(e) NA_character_))
                if (length(out) && !is.na(out[1])) return(parse_numeric_mb(out[1]))
                return(NA_real_)
        }
        
        if (grepl("windows", sys)) {
                # PowerShell: WorkingSet64（bytes）→ MB
                cmd <- sprintf('powershell -NoProfile -Command "(Get-Process -Id %d).WorkingSet64"', pid)
                out <- suppressWarnings(tryCatch(system(cmd, intern = TRUE), error=function(e) NA_character_))
                if (length(out) && !is.na(out[1])) {
                        bytes <- suppressWarnings(as.numeric(out[1]))
                        if (is.finite(bytes)) return(bytes / (1024^2))
                }
                # フォールバック
                val <- suppressWarnings(tryCatch(utils::memory.size(), error = function(e) NA_real_))
                if (!is.finite(val)) val <- NA_real_
                return(val)
        }
        
        NA_real_
}

# R(プロセス)の使用メモリ(MB)。pryr があれば使用、なければ OS 依存フォールバック
.mem_used_mb <- function() {
        if (requireNamespace("pryr", quietly = TRUE)) {
                val <- as.numeric(pryr::mem_used()) / (1024^2)
                if (!is.finite(val)) return(NA_real_)
                return(val)
        }
        if (.Platform$OS.type == "windows") {
                val <- suppressWarnings(tryCatch(utils::memory.size(), error = function(e) NA_real_))
                if (!is.finite(val)) val <- NA_real_
                return(val)
        }
        NA_real_
}

with_perf_log_extended <- function(tag, meta = list(),
                                   code,
                                   log_dir  = file.path(.get_base_dir(), "logs"),
                                   log_name = "perf_log_extended.csv") {
        code_sub <- substitute(code)
        
        # 前計測
        rss0 <- .get_rss_mb()
        mu0  <- .mem_used_mb()
        t0   <- proc.time()
        
        # 実行
        result <- eval(code_sub, envir = parent.frame())
        
        # 後計測
        pt   <- proc.time() - t0
        rss1 <- .get_rss_mb()
        mu1  <- .mem_used_mb()
        
        # 有限化
        nf <- function(x) { x <- suppressWarnings(as.numeric(x)); ifelse(is.finite(x), x, NA_real_) }
        
        row <- data.frame(
                timestamp   = as.character(Sys.time()),
                tag         = as.character(tag),         # "PE" or "IE"
                user_sec    = nf(pt[["user.self"]]),
                system_sec  = nf(pt[["sys.self"]]),
                elapsed_sec = nf(pt[["elapsed"]]),
                rss_mb_before       = nf(rss0),
                rss_mb_after        = nf(rss1),
                rss_mb_delta        = if (is.finite(rss0) && is.finite(rss1)) rss1 - rss0 else NA_real_,
                mem_used_mb_before  = nf(mu0),
                mem_used_mb_after   = nf(mu1),
                mem_used_mb_delta   = if (is.finite(mu0) && is.finite(mu1)) mu1 - mu0 else NA_real_,
                stringsAsFactors = FALSE
        )
        
        if (length(meta)) {
                meta_df <- as.data.frame(meta, optional = TRUE, stringsAsFactors = FALSE)
                row <- cbind(meta_df[rep(1, nrow(row)), , drop = FALSE], row)
        }
        
        dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
        log_path <- file.path(log_dir, log_name)
        
        if (!file.exists(log_path)) {
                utils::write.table(row, file = log_path, sep = ",", quote = TRUE,
                                   row.names = FALSE, col.names = TRUE, append = FALSE)
        } else {
                utils::write.table(row, file = log_path, sep = ",", quote = TRUE,
                                   row.names = FALSE, col.names = FALSE, append = TRUE)
        }
        
        invisible(result)
}

# 計測して CSV に1行追記しつつ、式を実行する
with_perf_log <- function(tag, meta = list(), code, log_dir = file.path(.get_base_dir(), "logs"), log_name = "perf_log.csv") {
        code_sub <- substitute(code)
        
        t0   <- proc.time()
        rss0 <- .get_rss_mb()
        mu0  <- .mem_used_mb()
        
        # 式を実行
        result <- eval(code_sub, envir = parent.frame())
        
        # 終了時点のメトリクス
        pt    <- proc.time() - t0
        rss1  <- .get_rss_mb()
        mu1   <- .mem_used_mb()
        
        row <- data.frame(
                timestamp   = as.character(Sys.time()),
                tag         = as.character(tag),
                user_sec    = unname(pt[["user.self"]]),
                system_sec  = unname(pt[["sys.self"]]),
                elapsed_sec = unname(pt[["elapsed"]]),
                rss_mb_before   = rss0,
                rss_mb_after    = rss1,
                rss_mb_delta    = if (is.na(rss0) || is.na(rss1)) NA_real_ else rss1 - rss0,
                mem_used_mb_before = mu0,
                mem_used_mb_after  = mu1,
                mem_used_mb_delta  = if (is.na(mu0) || is.na(mu1)) NA_real_ else mu1 - mu0,
                stringsAsFactors = FALSE
        )
        
        # メタ情報（DS, J, FileName, BSMethod, dt, thresholdMode, thresholdName など）を先頭に付加
        if (length(meta)) {
                meta_df <- as.data.frame(meta, optional = TRUE, stringsAsFactors = FALSE)
                # 列方向で結合（メタ列 → 計測列）
                row <- cbind(meta_df[rep(1, nrow(row)), , drop = FALSE], row)
        }
        
        # 追記
        dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
        log_path <- file.path(log_dir, log_name)
        utils::write.table(row, file = log_path, sep = ",", quote = TRUE,
                           row.names = FALSE, col.names = !file.exists(log_path), append = TRUE)
        
        invisible(result)
}

# === PE ===
PE <- function(DS, FileName, BSMethod, dt, thresholdMode, thresholdName, i) {
        base_dir <- .get_base_dir()
        BEout_Path <- file.path(base_dir, "BS-output")
        PEin_Path  <- file.path(BEout_Path, sprintf("%s/DS%s-BSCount-%s-%s-%s-%s(J=%s).dat",
                                                    FileName, DS, BSMethod, dt, thresholdMode, thresholdName, i))
        IEinpath   <- file.path(base_dir, "PE-output", FileName)
        .ensure_dir(IEinpath)
        dat_Path   <- create_pedat(DS, IEinpath, BSMethod, dt, thresholdMode, thresholdName, i)
        
        BS <- read.table(PEin_Path)
        bsNum <- nrow(BS)
        
        LSTu  <- vector("list", bsNum)
        LSTuC <- vector("list", bsNum)
        
        for (bsIndex in seq_len(bsNum)) {
                v <- as.numeric(BS[bsIndex, ])
                uPlot <- if (dt == "TI") {
                  Tipsh(v, thresholdMode, 1, i)
                } else if (dt == "H") {
                  Wse(v, "none", thresholdName, thresholdMode, 1, i)
                } else {
                  Wse(v, dt, thresholdName, thresholdMode, 1, i)
                }
                u <- uPlot$EstimationData
                # 丸めポリシーを統一（例：6桁）
                u <- round(u, 6)
                LSTu[[bsIndex]] <- u
                LSTuC[[bsIndex]] <- cumsum(u)
        }
        
        write.table(LSTu,  dat_Path$intensity_data,  append=FALSE, quote=FALSE, col.names=FALSE, row.names=FALSE)
        write.table(LSTuC, dat_Path$cumulative_data, append=FALSE, quote=FALSE, col.names=FALSE, row.names=FALSE)
}

# === IE ===
IE <- function(DS, FileName, BSMethod, dt, thresholdMode, thresholdName, i)
{
        if (!"package:e1071" %in% search()) suppressPackageStartupMessages(require(e1071))
        
        base_dir <- .get_base_dir()
        
        PE_data_Path    <- file.path(base_dir, "PE-output", FileName,
                                     sprintf("DS%s-PE-%s-k%s-%s-%s(J=%s).dat", DS, BSMethod, dt, thresholdMode, thresholdName, i))
        PECum_data_Path <- file.path(base_dir, "PE-output", FileName,
                                     sprintf("DS%s-PE-%s-r%s-%s-%s(J=%s).dat", DS, BSMethod, dt, thresholdMode, thresholdName, i))
        out_dir <- file.path(base_dir, "IE-output", FileName)
        .ensure_dir(out_dir)
        
        ds    <- as.data.frame(t(read.table(PE_data_Path)))
        dsCum <- as.data.frame(t(read.table(PECum_data_Path)))
        
        Num   <- ncol(ds)
        probs <- c(0.025, 0.975)
        pcnts <- c(0.25, 0.5, 0.75, 1.0)
        
        dsEachPnt    <- apply(ds,    2, sort, decreasing = FALSE)
        dsCumEachPnt <- apply(dsCum, 2, sort, decreasing = FALSE)
        
        final_path <- create_final(DS, i, out_dir, BSMethod, dt, thresholdMode, thresholdName)
        write.csv(t(dsEachPnt),    final_path$Pnt,    quote=FALSE, row.names = FALSE)
        write.csv(t(dsCumEachPnt), final_path$CumPnt, quote=FALSE, row.names = FALSE)
        
        MomIE  <- matrix(NA_real_, nrow = length(pcnts), ncol = 7)
        MomIEc <- matrix(NA_real_, nrow = length(pcnts), ncol = 7)
        
        # 四分位点の幅（UL-LL）
        w_ie_q  <- numeric(length(pcnts))
        w_iec_q <- numeric(length(pcnts))
        
        # 全時点の幅（UL-LL）を貯めるベクトル
        w_ie_all  <- numeric(Num)
        w_iec_all <- numeric(Num)
        
        # 0 の扱い：厳密ゼロを除外する例
        is_nz <- function(x, tol = 1e-12) abs(x) > tol
        
        .safe_moments <- function(x, LL, UL) {
                if (length(x) < 3L || all(is.na(x))) {
                        return(c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                                 round(LL,4), round(UL,4)))
                }
                c(
                        round(mean(x), 4),
                        round(median(x), 4),
                        round(var(x), 4),
                        round(e1071::skewness(x), 4),
                        round(e1071::kurtosis(x) + 3, 4),  # Pearson（正規=3）デフォルトだとtype=3,超過尖度
                        round(LL, 4),
                        round(UL, 4)
                )
        }
        
        # --- 四分位点での集計 ---
        for (k in seq_along(pcnts)) {
                point <- max(1, floor(Num * pcnts[k]))
                v  <- dsEachPnt[, point]
                vc <- dsCumEachPnt[, point]
                
                v_nz  <- v[is_nz(v)]
                vc_nz <- vc[is_nz(vc)]
                
                qv  <- if (length(v_nz)  > 0L) quantile(v_nz,  probs=probs, type=1, names=FALSE) else c(NA_real_, NA_real_)
                qvc <- if (length(vc_nz) > 0L) quantile(vc_nz, probs=probs, type=1, names=FALSE) else c(NA_real_, NA_real_)
                
                LL <- qv[1];  UL <- qv[2]
                LLc <- qvc[1]; ULc <- qvc[2]
                
                MomIE[k, ]  <- .safe_moments(as.numeric(v_nz),  LL,  UL)
                MomIEc[k, ] <- .safe_moments(as.numeric(vc_nz), LLc, ULc)
                
                w_ie_q[k]  <- as.numeric(UL  - LL)
                w_iec_q[k] <- as.numeric(ULc - LLc)
        }
        
        # --- 全時点での幅（UL-LL）を計算 ---
        for (j in seq_len(Num)) {
                vj  <- dsEachPnt[, j]
                vcj <- dsCumEachPnt[, j]
                vj  <- vj[is_nz(vj)]
                vcj <- vcj[is_nz(vcj)]
                
                qvj  <- if (length(vj)  > 0L) quantile(vj,  probs=probs, type=1, names=FALSE) else c(NA_real_, NA_real_)
                qvcj <- if (length(vcj) > 0L) quantile(vcj, probs=probs, type=1, names=FALSE) else c(NA_real_, NA_real_)
                
                w_ie_all[j]  <- as.numeric(qvj[2]  - qvj[1])
                w_iec_all[j] <- as.numeric(qvcj[2] - qvcj[1])
        }
        
        # 既存の7列出力はそのまま
        write.table(MomIE,  final_path$MomIE,  append=FALSE, quote=FALSE, col.names=FALSE, row.names=FALSE)
        write.table(MomIEc, final_path$MomIEc, append=FALSE, quote=FALSE, col.names=FALSE, row.names=FALSE)
        
        # 四分位点の幅（表 & 要約）
        ie_width_tbl  <- cbind(TimeFrac = pcnts, CI_Width = round(w_ie_q,  6))
        iec_width_tbl <- cbind(TimeFrac = pcnts, CI_Width = round(w_iec_q, 6))
        write.table(ie_width_tbl,  final_path$IEw,  append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
        write.table(iec_width_tbl, final_path$IEcw, append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
        
        ie_w_mean <- round(mean(w_ie_q,  na.rm=TRUE), 6)
        ie_w_var  <- round(var(w_ie_q,   na.rm=TRUE), 6)
        ie_w_sd   <- round(sd(w_ie_q,    na.rm=TRUE), 6)
        iec_w_mean <- round(mean(w_iec_q, na.rm=TRUE), 6)
        iec_w_var  <- round(var(w_iec_q,  na.rm=TRUE), 6)
        iec_w_sd   <- round(sd(w_iec_q,   na.rm=TRUE), 6)
        
        write.table(matrix(c(ie_w_mean,  ie_w_var,  ie_w_sd),  nrow=1),
                    final_path$IEw_sum,  append=FALSE, quote=FALSE, col.names=FALSE, row.names=FALSE)
        write.table(matrix(c(iec_w_mean, iec_w_var, iec_w_sd), nrow=1),
                    final_path$IEcw_sum, append=FALSE, quote=FALSE, col.names=FALSE, row.names=FALSE)
        
        # 全時点の幅の要約（mean, var, sd）を別ファイルに保存
        ie_all_mean <- round(mean(w_ie_all,  na.rm=TRUE), 6)
        ie_all_var  <- round(var(w_ie_all,   na.rm=TRUE), 6)
        ie_all_sd   <- round(sd(w_ie_all,    na.rm=TRUE), 6)
        iec_all_mean <- round(mean(w_iec_all, na.rm=TRUE), 6)
        iec_all_var  <- round(var(w_iec_all,  na.rm=TRUE), 6)
        iec_all_sd   <- round(sd(w_iec_all,   na.rm=TRUE), 6)
        
        write.table(matrix(c(ie_all_mean,  ie_all_var,  ie_all_sd),  nrow=1),
                    final_path$IEw_all_sum,  append=FALSE, quote=FALSE, col.names=FALSE, row.names=FALSE)
        write.table(matrix(c(iec_all_mean, iec_all_var, iec_all_sd), nrow=1),
                    final_path$IEcw_all_sum, append=FALSE, quote=FALSE, col.names=FALSE, row.names=FALSE)
}

# =========================================================
# 共通ヘルパ（パス統一・ディレクトリ作成）
# =========================================================
.get_base_dir <- function() {
        # 必要なら here::here() に置き換え可
        getwd()
}
.ensure_dir <- function(path) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# =========================================================
# creating file format (PE 出力ファイル名の生成)
# - 出力先ディレクトリが未作成でも安全に作成
# - 絶対パスで返す
# =========================================================
create_pedat <- function(DS, IEinpath, BSmethod, dt, thresholdMode, thresholdName, i)
{
        # 受け取った IEinpath を尊重（PE 内で base_dir/PE-output/FileName を渡す想定）
        .ensure_dir(IEinpath)
        
        dat_intensity  <- sprintf("/DS%s-PE-%s-k%s-%s-%s(J=%s).dat",
                                  DS, BSmethod, dt, thresholdMode, thresholdName, i)
        dat_cumulative <- sprintf("/DS%s-PE-%s-r%s-%s-%s(J=%s).dat",
                                  DS, BSmethod, dt, thresholdMode, thresholdName, i)
        
        intensity_data  <- paste0(IEinpath, dat_intensity)
        cumulative_data <- paste0(IEinpath, dat_cumulative)
        
        list(intensity_data = intensity_data, cumulative_data = cumulative_data)
}

# =========================================================
# creating file format (IE 出力ファイル名の生成)
# - 出力先ディレクトリが未作成でも安全に作成
# - 絶対パスで返す
# =========================================================
create_final <- function(DS,i, directory_path,BSMethod,dt,thresholdMode,thresholdName)
{
        .ensure_dir(directory_path)
        
        file_name_Pnt    <- sprintf("/DS%s-IEEachPnt-%s-k%s-%s-%s(J=%s).csv",
                                    DS, BSMethod, dt, thresholdMode, thresholdName, i)
        file_name_CumPnt <- sprintf("/DS%s-IEEachPnt-%s-r%s-%s-%s(J=%s).csv",
                                    DS, BSMethod, dt, thresholdMode, thresholdName, i)
        file_name_MomIE  <- sprintf("/DS%s-IE-%s-%s-%s-%s(J=%s).dat",
                                    DS, BSMethod, dt, thresholdMode, thresholdName, i)
        file_name_MomIEc <- sprintf("/DS%s-IEc-%s-%s-%s-%s(J=%s).dat",
                                    DS, BSMethod, dt, thresholdMode, thresholdName, i)
        
        # 四分位点の幅 & 要約
        file_name_IEw      <- sprintf("/DS%s-IEwidth-%s-%s-%s-%s(J=%s).dat",
                                      DS, BSMethod, dt, thresholdMode, thresholdName, i)
        file_name_IEcw     <- sprintf("/DS%s-IEcwidth-%s-%s-%s-%s(J=%s).dat",
                                      DS, BSMethod, dt, thresholdMode, thresholdName, i)
        file_name_IEw_sum  <- sprintf("/DS%s-IEwidth-summary-%s-%s-%s-%s(J=%s).dat",
                                      DS, BSMethod, dt, thresholdMode, thresholdName, i)
        file_name_IEcw_sum <- sprintf("/DS%s-IEcwidth-summary-%s-%s-%s-%s(J=%s).dat",
                                      DS, BSMethod, dt, thresholdMode, thresholdName, i)
        
        # 全時点の幅の要約（mean, var, sd）
        file_name_IEw_all_sum  <- sprintf("/DS%s-IEwidth-all-summary-%s-%s-%s-%s(J=%s).dat",
                                          DS, BSMethod, dt, thresholdMode, thresholdName, i)
        file_name_IEcw_all_sum <- sprintf("/DS%s-IEcwidth-all-summary-%s-%s-%s-%s(J=%s).dat",
                                          DS, BSMethod, dt, thresholdMode, thresholdName, i)
        
        Pnt    <- paste0(directory_path, file_name_Pnt)
        CumPnt <- paste0(directory_path, file_name_CumPnt)
        MomIE  <- paste0(directory_path, file_name_MomIE)
        MomIEc <- paste0(directory_path, file_name_MomIEc)
        IEw      <- paste0(directory_path, file_name_IEw)
        IEcw     <- paste0(directory_path, file_name_IEcw)
        IEw_sum  <- paste0(directory_path, file_name_IEw_sum)
        IEcw_sum <- paste0(directory_path, file_name_IEcw_sum)
        IEw_all_sum  <- paste0(directory_path, file_name_IEw_all_sum)
        IEcw_all_sum <- paste0(directory_path, file_name_IEcw_all_sum)
        
        list(Pnt=Pnt, CumPnt=CumPnt, MomIE=MomIE, MomIEc=MomIEc,
             IEw=IEw, IEcw=IEcw, IEw_sum=IEw_sum, IEcw_sum=IEcw_sum,
             IEw_all_sum=IEw_all_sum, IEcw_all_sum=IEcw_all_sum)
}

# =========================================================
# CI（後方互換 + 改善） 未使用（互換）
# - 推奨：順序統計量（Non0）から直接分位点を返す（経験分布の定義通り）
# - 旧：ヒストグラム累積を渡していた場合にも、できるだけ安全に動作
# =========================================================
CI <- function(Distri, h, Non0, x) {
        # 新方式：Non0 が数列ベクトルなら、その経験分布から分位点を返す
        if (!missing(Non0) && is.numeric(Non0) && length(Non0) >= 1) {
                return(as.numeric(quantile(as.numeric(Non0), probs = h, type = 1, names = FALSE)))
        }
        
        # 旧方式へのフォールバック（非推奨）
        # Distri: 累積（0..1）を想定、x: ビン幅、Non0[1] を下端としていたが曖昧
        warning("CI(): ヒストグラム由来の境界推定は非推奨です。可能なら Non0(サンプル) を渡してください。")
        dis <- Distri
        com <- h
        idx_lo <- max(1L, sum(dis <= com))
        idx_hi <- min(length(dis), idx_lo + 1L)
        rs <- ifelse(idx_lo == 1L, dis[1L],
                     ifelse(abs(com - dis[idx_lo]) <= abs(dis[idx_hi] - com), dis[idx_lo], dis[idx_hi]))
        # 下端推定がないため、元の仕様通りに x * 位置 + 0 を返す（下端0仮定）
        # 実質的に数量スケールの復元ができないため、非推奨扱い
        which_idx <- match(rs, dis)
        as.numeric(which_idx) * (ifelse(missing(x), 1, x))
}

# =========================================================
# BSintervals（ラッパ）
# - パス統一の上で PE → IE を呼び出し
# - 主要出力パスを返して後続処理にも使いやすく
# =========================================================
BSintervals <- function(DS, FileName, BSMethod, dt, thresholdMode, thresholdName, i){
        base_dir <- .get_base_dir()
        
        meta <- list(
                dataset = DS, J = i, FileName = FileName,
                BSMethod = BSMethod, dt = dt,
                thresholdMode = thresholdMode, thresholdName = thresholdName
        )
        
        with_perf_log_extended("PE", meta, {
                PE(DS, FileName, BSMethod, dt, thresholdMode, thresholdName, i)
        })
        
        with_perf_log_extended("IE", meta, {
                IE(DS, FileName, BSMethod, dt, thresholdMode, thresholdName, i)
        })
        
        list(
                pe_dir = file.path(base_dir, "PE-output", FileName),
                ie_dir = file.path(base_dir, "IE-output", FileName)
        )
}
