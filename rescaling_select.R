# =========================================================
# メイン1（点推定）までの NHPP 診断スクリプト
# 出力: CSV（要約）, 図（ggplot2; uのヒスト/ECDF/QQ, −log10 pバー）
# 保存先: output/_summaries/<GroupDir>/<FileName>/
# =========================================================

# 必要パッケージ
suppressPackageStartupMessages({
        library(ggplot2)
})
suppressPackageStartupMessages({ library(ggplot2) })

.ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

# ---------- 読み込みヘルパ ----------
read_counts_vec <- function(path_counts) {
        # DS/DS1.txt のように (bin, count) 2列（空白/タブ区切り）を想定
        df <- utils::read.table(path_counts, header = FALSE, sep = "",
                                check.names = FALSE, comment.char = "",
                                quote = "", stringsAsFactors = FALSE)
        as.numeric(df[[2]])
}

read_lambda_vec <- function(path_lambda) {
        # output/.../DS1-kA1-h-ut(J=2).txt のように 1列CSV（先頭ヘッダあり）
        df <- utils::read.table(path_lambda, header = TRUE, sep = ",",
                                check.names = FALSE, comment.char = "",
                                stringsAsFactors = FALSE)
        as.numeric(df[[1]])
}

# ---------- NHPP 診断（メイン1のみ） ----------
nhpp_diagnostics_main1 <- function(path_counts, path_lambda, dt = 1, min_exp = 5, jitter_eps = 1e-12) {
        counts <- read_counts_vec(path_counts)
        lam    <- read_lambda_vec(path_lambda)
        
        if (length(counts) != length(lam)) {
                stop(sprintf("length mismatch: counts=%d vs lambda=%d", length(counts), length(lam)))
        }
        
        # 1) χ²：期待= λ * Δt、期待が小さいビンは右へ順次マージ
        exp_counts <- lam * dt
        O <- c(); E <- c(); accO <- 0; accE <- 0
        for (j in seq_along(counts)) {
                accO <- accO + counts[j]; accE <- accE + exp_counts[j]
                if (accE >= min_exp || j == length(counts)) {
                        O <- c(O, accO); E <- c(E, accE); accO <- 0; accE <- 0
                }
        }
        df_chi <- length(O) - 1L
        chisq  <- sum((O - E)^2 / pmax(E, .Machine$double.eps))
        p_chi  <- stats::pchisq(chisq, df = df_chi, lower.tail = FALSE)
        
        # 2) 近似 time-rescaling → U(0,1)
        u <- 1 - exp(-(lam * dt))
        # ties回避の微小ジッタ
        u <- pmin(pmax(u + stats::runif(length(u), -jitter_eps, jitter_eps), 0), 1)
        u <- u[is.finite(u) & u > 0 & u < 1]
        
        ks_res <- stats::ks.test(u, "punif")
        
        # AD / CvM（goftest があれば）
        if (!requireNamespace("goftest", quietly = TRUE)) {
                ad_p  <- NA_real_; cvm_p <- NA_real_
                ad_msg <- "goftest not installed; AD skipped"
        } else {
                ad_p  <- tryCatch(goftest::ad.test(u, null = "punif")$p.value,  error = function(e) NA_real_)
                cvm_p <- tryCatch(goftest::cvm.test(u, null = "punif")$p.value, error = function(e) NA_real_)
                ad_msg <- NULL
        }
        
        list(
                counts = counts,
                lambda = lam,
                u = u,
                chisq = list(stat = chisq, df = df_chi, p.value = p_chi, observed = O, expected = E),
                ks_p  = ks_res$p.value,
                ad_p  = ad_p,
                cvm_p = cvm_p,
                note  = ad_msg
        )
}

# ---------- 図とCSVの保存（白背景固定版） ----------
save_main1_outputs <- function(diag, out_dir, tag) {
        .ensure_dir(out_dir)
        
        # 要約テーブルをCSVで
        summary_df <- data.frame(
                tag = tag,
                chisq_stat = diag$chisq$stat,
                chisq_df   = diag$chisq$df,
                chisq_p    = diag$chisq$p.value,
                ks_p       = diag$ks_p,
                ad_p       = diag$ad_p,
                cvm_p      = diag$cvm_p,
                stringsAsFactors = FALSE
        )
        utils::write.csv(summary_df, file.path(out_dir, paste0("main1_nhpp_summary_", tag, ".csv")), row.names = FALSE)
        
        # 共通テーマ（白背景）
        white_bg <- theme_minimal(base_size = 14) +
                theme(
                        panel.background = element_rect(fill = "white", colour = NA),
                        plot.background  = element_rect(fill = "white", colour = NA)
                )
        
        # 1) u のヒストグラム
        if (length(diag$u)) {
                p_hist <- ggplot(data.frame(u = diag$u), aes(x = u)) +
                        geom_histogram(bins = max(10, floor(sqrt(length(diag$u)))), boundary = 0, closed = "right",
                                       color = "black", fill = "steelblue", alpha = 0.8) +
                        geom_hline(yintercept = length(diag$u)/10, linetype = 2, linewidth = 0.4, alpha = 0.7) +
                        labs(title = sprintf("Histogram of u (rescaling) — %s", tag), x = "u ~ U(0,1)", y = "count") +
                        white_bg
                ggsave(file.path(out_dir, paste0("main1_u_hist_", tag, ".png")),
                       p_hist, width = 6, height = 4, dpi = 120, bg = "white")
        }
        
        # 2) u の ECDF vs Uniform の45°線
        if (length(diag$u)) {
                df_ecdf <- data.frame(u = sort(diag$u), Fhat = ecdf(diag$u)(sort(diag$u)))
                p_ecdf <- ggplot(df_ecdf, aes(x = u, y = Fhat)) +
                        geom_step() +
                        geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4) +
                        labs(title = sprintf("ECDF of u vs Uniform — %s", tag), x = "u", y = "F̂(u)") +
                        white_bg
                ggsave(file.path(out_dir, paste0("main1_u_ecdf_", tag, ".png")),
                       p_ecdf, width = 6, height = 4, dpi = 120, bg = "white")
        }
        
        # 3) u の QQ（Uniform(0,1)）
        if (length(diag$u)) {
                n <- length(diag$u)
                u_sorted <- sort(diag$u)
                q_theory <- (1:n)/(n+1)  # 理論分位
                dfqq <- data.frame(theory = q_theory, sample = u_sorted)
                p_qq <- ggplot(dfqq, aes(x = theory, y = sample)) +
                        geom_point(size = 1.4, alpha = 0.8) +
                        geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4) +
                        labs(title = sprintf("QQ plot: Uniform vs u — %s", tag), x = "Uniform quantile", y = "Sample quantile") +
                        white_bg
                ggsave(file.path(out_dir, paste0("main1_u_qq_", tag, ".png")),
                       p_qq, width = 6, height = 4, dpi = 120, bg = "white")
        }
        
        # 4) −log10 p（KS/AD/CvM）バー
        pvals <- c(KS = diag$ks_p, AD = diag$ad_p, CvM = diag$cvm_p)
        logp  <- -log10(pmax(pvals, .Machine$double.xmin))
        dfbar <- data.frame(test = names(pvals), logp = as.numeric(logp))
        p_bar <- ggplot(dfbar, aes(x = test, y = logp, fill = test)) +
                geom_col() +
                geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.4) +
                labs(title = sprintf("-log10(p): KS / AD / CvM — %s", tag), x = "test", y = "-log10(p)") +
                white_bg +
                theme(legend.position = "none")
        ggsave(file.path(out_dir, paste0("main1_logp_bar_", tag, ".png")),
               p_bar, width = 6, height = 4, dpi = 120, bg = "white")
        
        # テキストログ（オプション）
        sink(file.path(out_dir, paste0("main1_model_level_", tag, ".txt")))
        cat("== NHPP diagnostics (Main1 only) ==\n\n")
        print(diag$chisq)
        cat("\nKS p:", diag$ks_p, "\nAD p:", diag$ad_p, "\nCvM p:", diag$cvm_p, "\n")
        if (!is.null(diag$note)) cat("\nNote:", diag$note, "\n")
        sink(NULL)
        
        invisible(summary_df)
}

# ---------- 1ケース実行 ----------
run_main1_diagnostics <- function(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i,
                                  out_root = "output/_summaries", dt = 1, seed = 123) {
        set.seed(seed)
        tag <- sprintf("DS%s_%s_%s_i%s_k%s", DS, thresholdMode, thresholdName, i, dt_code)
        
        # ★ 出力ディレクトリ：FileName が空なら階層を落とす
        out_dir <- if (!is.null(FileName) && nzchar(FileName)) {
                file.path(out_root, GroupDir, FileName)
        } else {
                file.path(out_root, GroupDir)
        }
        .ensure_dir(out_dir)
        
        path_counts <- file.path("DS", sprintf("DS%s.txt", DS))
        
        # ★ λファイルの場所：FileName が空なら省略
        lambda_dir <- if (!is.null(FileName) && nzchar(FileName)) {
                file.path("output", GroupDir, FileName)
        } else {
                file.path("output", GroupDir)
        }
        path_lambda <- file.path(
                lambda_dir,
                sprintf("DS%s-k%s-%s-%s(J=%s).txt", DS, dt_code, thresholdMode, thresholdName, i)
        )
        
        if (!file.exists(path_counts)) stop("counts file not found: ", path_counts)
        if (!file.exists(path_lambda)) stop("lambda (Main1) file not found: ", path_lambda)
        
        diag <- nhpp_diagnostics_main1(path_counts, path_lambda, dt = dt)
        save_main1_outputs(diag, out_dir, tag)
        invisible(diag)
}


# ---------- 複数ケース一括 ----------
batch_main1_diagnostics <- function(cases,
                                    out_root = "output/_summaries",
                                    dt = 1, seed = 123) {
        set.seed(seed)
        all_rows <- list()
        for (r in seq_len(nrow(cases))) {
                DS <- as.character(cases$DS[r])
                GroupDir <- as.character(cases$GroupDir[r])
                FileName <- as.character(cases$FileName[r])
                dt_code <- as.character(cases$dt_code[r])
                thresholdMode <- as.character(cases$thresholdMode[r])
                thresholdName <- as.character(cases$thresholdName[r])
                i <- as.character(cases$i[r])
                
                message(sprintf("Main1 diagnostics: DS=%s, %s/%s, k=%s, %s/%s, i=%s",
                                DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i))
                
                # 1ケースずつ実行（保存も中で行う）
                summ <- tryCatch({
                        diag <- run_main1_diagnostics(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i,
                                                      out_root = out_root, dt = dt, seed = seed)
                        # 戻り値の summary から1行化
                        data.frame(
                                DS = DS, GroupDir = GroupDir, FileName = FileName,
                                dt_code = dt_code, thresholdMode = thresholdMode, thresholdName = thresholdName, i = i,
                                chisq_stat = diag$chisq$stat, chisq_df = diag$chisq$df, chisq_p = diag$chisq$p.value,
                                ks_p = diag$ks_p, ad_p = diag$ad_p, cvm_p = diag$cvm_p,
                                stringsAsFactors = FALSE
                        )
                }, error = function(e) {
                        warning("Main1 diagnostics failed for DS=", DS, ", i=", i, " : ", conditionMessage(e))
                        data.frame(
                                DS = DS, GroupDir = GroupDir, FileName = FileName,
                                dt_code = dt_code, thresholdMode = thresholdMode, thresholdName = thresholdName, i = i,
                                chisq_stat = NA_real_, chisq_df = NA_real_, chisq_p = NA_real_,
                                ks_p = NA_real_, ad_p = NA_real_, cvm_p = NA_real_,
                                stringsAsFactors = FALSE
                        )
                })
                all_rows[[length(all_rows) + 1L]] <- summ
        }
        
        big <- do.call(rbind, all_rows)
        # DS/FileName ごとに総括CSVを保存
        # 既存: combos <- unique(big[, c("GroupDir","FileName")])
        combos <- unique(big[, c("GroupDir","FileName")])
        
        for (k in seq_len(nrow(combos))) {
                gd <- combos$GroupDir[k]
                fn <- combos$FileName[k]
                
                # ★ FileName が空なら階層を落とす
                out_dir <- if (!is.null(fn) && nzchar(fn)) {
                        file.path(out_root, gd, fn)
                } else {
                        file.path(out_root, gd)
                }
                .ensure_dir(out_dir)
                
                sel <- if (!is.null(fn) && nzchar(fn)) {
                        big$GroupDir == gd & big$FileName == fn
                } else {
                        big$GroupDir == gd & (is.na(big$FileName) | big$FileName == "")
                }
                
                utils::write.csv(
                        big[sel, ],
                        file.path(out_dir, "main1_all_cases_summary.csv"),
                        row.names = FALSE
                )
        }
        invisible(big)
}

make_main1_heatmaps_focus <- function(csv_path,
                                      targets,
                                      out_root = NULL,
                                      na_color = "grey90",
                                      logp_cap = 50) {
        if (!file.exists(csv_path)) stop("CSV not found: ", csv_path)
        
        df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
        
        # 必要列チェック
        req <- c("DS","GroupDir","FileName","i","thresholdMode","thresholdName","chisq_p","ks_p","ad_p","cvm_p")
        missing <- setdiff(req, names(df))
        if (length(missing)) stop("CSV missing columns: ", paste(missing, collapse=", "))
        
        # 絞り込み
        df <- subset(df,
                     (FileName=="A3" & thresholdMode=="h" & thresholdName=="ut") |
                             (FileName=="B2" & thresholdMode=="h" & thresholdName=="ut") |
                             (FileName==""   & GroupDir=="NDT_WSE" & thresholdMode=="h" & thresholdName=="ldt"))
        
        if (!nrow(df)) stop("対象データがありません")
        
        # 長データ化
        to_long <- function(d) {
                keep <- d[, c("DS","i","GroupDir","FileName","thresholdMode","thresholdName","ks_p","ad_p","cvm_p")]
                long <- reshape(
                        data = keep,
                        varying = list(7:9), v.names = "pval",
                        timevar = "test", times = c("KS","AD","CvM"),
                        direction = "long"
                )
                long$logp <- -log10(pmax(long$pval, .Machine$double.xmin))
                long$logp <- pmin(long$logp, logp_cap)  # 上限カット
                long$i_num <- suppressWarnings(as.numeric(as.character(long$i)))
                long$test  <- factor(long$test, levels = c("KS","AD","CvM"))
                long
        }
        
        df_long <- to_long(df)
        if (!nrow(df_long)) stop("データがありません")
        
        if (is.null(out_root)) out_root <- dirname(csv_path)
        .ensure_dir(out_root)
        
        # GroupDir+FileName ごとに出力
        combos <- unique(df_long[,c("GroupDir","FileName","thresholdMode","thresholdName")])
        for (k in seq_len(nrow(combos))) {
                sub <- subset(df_long,
                              GroupDir==combos$GroupDir[k] &
                                      FileName==combos$FileName[k] &
                                      thresholdMode==combos$thresholdMode[k] &
                                      thresholdName==combos$thresholdName[k])
                
                if (!nrow(sub)) next
                
                fn <- if (nzchar(combos$FileName[k])) combos$FileName[k] else "NDT"
                title_txt <- sprintf("Main1 NHPP diagnostics — capped −log10(p) [%s/%s] (%s/%s)",
                                     combos$GroupDir[k], fn,
                                     combos$thresholdMode[k], combos$thresholdName[k])
                
                p <- ggplot(sub, aes(x = test, y = i_num, fill = logp)) +
                        geom_tile() +
                        scale_fill_gradient(low="white", high="purple", na.value=na_color,
                                            limits=c(0, logp_cap)) +
                        facet_wrap(~ DS, nrow=1, scales="free_y") +
                        labs(title=title_txt,
                             x="test (KS/AD/CvM)", y="i", fill="-log10(p) (capped)") +
                        theme_minimal(base_size=14) +
                        theme(panel.background = element_rect(fill="white", colour=NA),
                              plot.background  = element_rect(fill="white", colour=NA))
                
                ggsave(file.path(out_root,
                                 sprintf("main1_focus_heatmap_%s_%s_%s_%s.png",
                                         combos$GroupDir[k], fn,
                                         combos$thresholdMode[k], combos$thresholdName[k])),
                       p, width=8, height=5, dpi=150, bg="white")
        }
        
        invisible(TRUE)
}

# --- 3条件(A3/B2/NDT-H)に限定して、正しいCSVを読み→出力は output/_summaries ---
make_main1_heatmaps_focus3 <- function(
                ans_csv = file.path("output/_summaries","DT_Ans_WSE","A3","main1_all_cases_summary.csv"),
                bar_csv = file.path("output/_summaries","DT_Bar_WSE","B2","main1_all_cases_summary.csv"),
                ndt_csv = file.path("output/_summaries","NDT_WSE","main1_all_cases_summary.csv"),
                out_root = file.path("output","_summaries"),
                na_color = "grey90",
                logp_cap = 50
){
        .ensure_dir(out_root)
        
        # 1) 読み込み（存在するものだけ）
        read_ok <- function(p){
                if (!file.exists(p)) {
                        message("CSV not found (skip): ", p); return(NULL)
                }
                utils::read.csv(p, stringsAsFactors = FALSE)
        }
        df_list <- Filter(Negate(is.null), list(read_ok(ans_csv), read_ok(bar_csv), read_ok(ndt_csv)))
        if (!length(df_list)) stop("対象CSVが見つかりませんでした。")
        
        df <- do.call(rbind, df_list)
        df$FileName[is.na(df$FileName)] <- ""
        
        # 必須列チェック
        req <- c("DS","GroupDir","FileName","i","thresholdMode","thresholdName",
                 "chisq_p","ks_p","ad_p","cvm_p")
        miss <- setdiff(req, names(df))
        if (length(miss)) stop("CSV missing columns: ", paste(miss, collapse=", "))
        
        # 2) 3条件でフィルタ
        df_focus <- subset(
                df,
                (GroupDir=="DT_Ans_WSE" & FileName=="A3" & thresholdMode=="h" & thresholdName=="ut") |
                        (GroupDir=="DT_Bar_WSE" & FileName=="B2" & thresholdMode=="h" & thresholdName=="ut") |
                        (GroupDir=="NDT_WSE"    & (is.na(FileName) | FileName=="") &
                                 thresholdMode=="h" & thresholdName=="ldt")
        )
        if (!nrow(df_focus)) stop("3条件に一致する行がありません。")
        
        # 3) ロング化
        to_long <- function(d){
                keep <- d[, c("DS","i","GroupDir","FileName","thresholdMode","thresholdName","ks_p","ad_p","cvm_p")]
                long <- reshape(
                        data = keep,
                        varying = list(7:9), v.names = "pval",
                        timevar = "test", times = c("KS","AD","CvM"),
                        direction = "long"
                )
                long$logp  <- -log10(pmax(long$pval, .Machine$double.xmin))
                long$logp  <- pmin(long$logp, logp_cap)  # 上限カット
                long$i_num <- suppressWarnings(as.numeric(as.character(long$i)))
                long$test  <- factor(long$test, levels = c("KS","AD","CvM"))
                long
        }
        long <- to_long(df_focus)
        long <- long[is.finite(long$i_num), , drop = FALSE]
        
        # 4) GroupDir + FileName(空なら省略) + しきい値ごとに出力（出力先は out_root）
        combos <- unique(long[, c("GroupDir","FileName","thresholdMode","thresholdName")])
        
        for (k in seq_len(nrow(combos))) {
                gd <- combos$GroupDir[k]
                fn <- combos$FileName[k]           # "" が NDT
                tm <- combos$thresholdMode[k]
                tn <- combos$thresholdName[k]
                
                # long から必要条件でサブセット（★ しきい値も含めて絞る）
                sub <- subset(
                        long,
                        GroupDir == gd &
                                (FileName == fn | (fn == "" & FileName == "")) &
                                thresholdMode == tm & thresholdName == tn
                )
                if (!nrow(sub)) next
                
                # ラベル（NDT は FileName を省略表示）
                fn_label <- if (nzchar(fn)) fn else "NDT"
                
                title_txt <- sprintf(
                        "Main1 NHPP diagnostics — capped −log10(p) [%s/%s] (%s/%s)",
                        gd, fn_label, tm, tn
                )
                
                p <- ggplot(sub, aes(x = test, y = i_num, fill = logp)) +
                        geom_tile() +
                        scale_fill_gradient(low = "white", high = "purple", na.value = na_color,
                                            limits = c(0, logp_cap)) +
                        facet_wrap(~ DS, nrow = 1, scales = "free_y") +
                        labs(title = title_txt, x = "test (KS / AD / CvM)", y = "i",
                             fill = "-log10(p) (capped)") +
                        theme_minimal(base_size = 14) +
                        theme(panel.background = element_rect(fill = "white", colour = NA),
                              plot.background  = element_rect(fill = "white", colour = NA))
                
                ggsave(file.path(out_root,
                                 sprintf("main1_focus_heatmap_%s_%s_%s_%s.png",
                                         gd, fn_label, tm, tn)),
                       p, width = 8, height = 5, dpi = 150, bg = "white")
        }
        
        invisible(TRUE)
}

make_main1_heatmaps_focus3(
        ans_csv = "output/_summaries/DT_Ans_WSE/A3/main1_all_cases_summary.csv",
        bar_csv = "output/_summaries/DT_Bar_WSE/B2/main1_all_cases_summary.csv",
        ndt_csv = "output/_summaries/NDT_WSE/main1_all_cases_summary.csv",
        out_root = "output/_summaries"  # ←出力はここに固定
)