library(ggplot2)

.ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
.boot_method_order <- c("COS","GS-","GS+","OrdSta","ThiArr","ThiInt")

#############################
# パスビルダ
#############################
build_paths <- function(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i, BSMethod=NULL) {
        # FileName を安全に正規化（NA -> ""）
        FileName <- as.character(FileName)
        if (is.na(FileName)) FileName <- ""
        
        # FileName が空なら中間フォルダを挟まない
        out_dir <- if (nzchar(FileName)) file.path("output",   GroupDir, FileName) else file.path("output",   GroupDir)
        ie_dir  <- if (nzchar(FileName)) file.path("IE-output", GroupDir, FileName) else file.path("IE-output", GroupDir)
        
        list(
                counts = file.path("DS", sprintf("DS%s.txt", DS)),
                # メイン1の観測（点推定λ）パス
                obs    = file.path(out_dir, sprintf("DS%s-k%s-%s-%s(J=%s).txt",
                                                    DS, dt_code, thresholdMode, thresholdName, i)),
                # IEEachPnt：BSMethod が指定されている時だけ
                ie_each = if (!is.null(BSMethod)) {
                        file.path(ie_dir, sprintf("DS%s-IEEachPnt-%s-k%s-%s-%s(J=%s).csv",
                                                  DS, BSMethod, dt_code, thresholdMode, thresholdName, i))
                } else {
                        NULL
                }
        )
}

#############################
# 読み込みヘルパ（安全版）
#############################
read_counts_vec <- function(path_counts) {
        df <- utils::read.table(path_counts, header = FALSE, sep = "",
                                check.names = FALSE, comment.char = "",
                                quote = "", stringsAsFactors = FALSE)
        as.numeric(df[[2]])
}

read_obs_hat_vec <- function(path_obs) {
        df <- utils::read.table(path_obs, header = TRUE, sep = ",",
                                check.names = FALSE, comment.char = "",
                                stringsAsFactors = FALSE)
        as.numeric(df[[1]])
}

# IEEachPnt：行=時点、列=ブート標本
read_ie_matrix_rows_time <- function(path_ie_each_csv) {
        ie <- utils::read.csv(path_ie_each_csv, header = TRUE,
                              check.names = FALSE, stringsAsFactors = FALSE)
        # 全NA列除去
        is_all_na <- vapply(ie, function(col) all(is.na(col)), logical(1))
        if (any(is_all_na)) ie <- ie[!is_all_na]
        # 先頭列が 1..nrow の純インデックスなら落とす
        to_num <- function(x) suppressWarnings(as.numeric(x))
        if (ncol(ie) >= 1) {
                v1 <- to_num(ie[[1]])
                if (all(is.finite(v1)) && length(v1) == nrow(ie) &&
                    all(abs(v1 - seq_len(nrow(ie))) < 1e-9)) {
                        ie <- ie[-1]
                }
        }
        # 数値化できない列を除去し、数値行列へ
        num_ok <- vapply(ie, function(col) { y <- to_num(col); any(is.finite(y)) }, logical(1))
        ie <- ie[num_ok]
        as.matrix(data.frame(lapply(ie, to_num)))
}

#############################
# 検定ユーティリティ
#############################
test_nhpp_chisq <- function(path_counts, path_obs_lambda, dt = 1, min_exp = 5) {
        counts <- read_counts_vec(path_counts)
        lam    <- read_obs_hat_vec(path_obs_lambda)
        stopifnot(length(counts) == length(lam))
        exp_counts <- lam * dt
        O <- c(); E <- c(); accO <- 0; accE <- 0
        for (j in seq_along(counts)) {
                accO <- accO + counts[j]; accE <- accE + exp_counts[j]
                if (accE >= min_exp || j == length(counts)) {
                        O <- c(O, accO); E <- c(E, accE); accO <- 0; accE <- 0
                }
        }
        df <- length(O) - 1
        chisq <- sum((O - E)^2 / pmax(E, .Machine$double.eps))
        pval  <- stats::pchisq(chisq, df = df, lower.tail = FALSE)
        list(stat = chisq, df = df, p.value = pval, observed = O, expected = E)
}

test_rescaling_uniform <- function(path_obs_lambda, dt = 1, eps = 1e-12) {
        lam <- read_obs_hat_vec(path_obs_lambda)
        u <- 1 - exp(-(lam * dt))
        u <- pmin(pmax(u + stats::runif(length(u), -eps, eps), 0), 1)  # ジッタ
        u <- u[is.finite(u) & u > 0 & u < 1]
        stats::ks.test(u, "punif")
}

.rand_pit_one <- function(sample_vec, x) {
  s <- sample_vec[is.finite(sample_vec)]
  m <- length(s)
  if (m == 0L || !is.finite(x)) return(NA_real_)
  Fx_minus <- sum(s <  x) / m
  Fx       <- sum(s <= x) / m
  
  if (Fx_minus == Fx) {
    # 全サンプルがxに等しいケース → (0,1)からランダムに返す
    return(stats::runif(1L))
  } else {
    return(Fx_minus + stats::runif(1L) * (Fx - Fx_minus))
  }
}

compute_pit_vector <- function(path_ie_each_csv, path_obs_lambda, seed = 123) {
        obs <- read_obs_hat_vec(path_obs_lambda)           # 長さ=時点数
        ie  <- read_ie_matrix_rows_time(path_ie_each_csv)  # nrow=時点, ncol=標本
        L <- min(nrow(ie), length(obs))
        ie <- ie[seq_len(L), , drop=FALSE]; obs <- obs[seq_len(L)]
        set.seed(seed)
        pit <- vapply(seq_len(L), function(j) .rand_pit_one(ie[j, ], obs[j]), numeric(1))
        pit[is.finite(pit) & pit > 0 & pit < 1]
}

ks_pit_with_jitter <- function(pit, eps = 1e-12) {
  # 前処理：有限 & (0,1) に限定
  pit <- pit[is.finite(pit) & pit > 0 & pit < 1]
  
  # ★ ここが重要：サンプル不足なら NA を返す
  if (length(pit) < 2L) {
    out <- list(statistic = NA_real_, p.value = NA_real_, method = "KS (insufficient PIT)",
                data.name = "pit")
    class(out) <- "htest"
    return(out)
  }
  
  pit2 <- pmin(pmax(pit + stats::runif(length(pit), -eps, eps), 0), 1)
  stats::ks.test(pit2, "punif")
}

chisq_uniformity <- function(cov_vec, nbins = 10) {
  O <- table(cut(cov_vec, breaks = seq(0, 1, length.out = nbins + 1), include.lowest = TRUE))
  O <- as.numeric(O)
  
  # ★ ここを追加: 全ゼロなら NA を返す
  if (all(O == 0)) {
    out <- list(statistic = NA_real_, p.value = NA_real_, method = "Chi-sq (all zero O)", data.name = "coverage")
    class(out) <- "htest"
    return(out)
  }
  
  stats::chisq.test(O, p = rep(1/nbins, nbins), simulate.p.value = TRUE)
}

test_coverage_from_ieeach <- function(path_ie_each_csv, path_obs_lambda,
                                      probs = c(0.025, 0.975)) {
        ie  <- read_ie_matrix_rows_time(path_ie_each_csv)
        obs <- read_obs_hat_vec(path_obs_lambda)
        L <- min(nrow(ie), length(obs))
        ie <- ie[seq_len(L), , drop=FALSE]; obs <- obs[seq_len(L)]
        Q <- t(apply(ie, 1, function(row) {
                row <- row[is.finite(row)]
                if (!length(row)) return(c(NA_real_, NA_real_))
                stats::quantile(row, probs = probs, type = 1, names = FALSE)
        }))
        LL <- Q[,1]; UL <- Q[,2]
        hit <- as.integer(LL <= obs & obs <= UL)
        k <- sum(hit, na.rm = TRUE)
        n <- sum(is.finite(LL) & is.finite(UL) & is.finite(obs))
        nominal <- 1 - 2 * probs[1]
        list(k = k, n = n, coverage = if (n > 0) k / n else NA_real_,
             binom_p = stats::binom.test(k, n, p = nominal)$p.value)
}

#############################
# メイン：複数BSMethodの比較
#############################
compare_boot_methods <- function(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i,
                                 methods = c("COS","GS-","GS+","OrdSta","ThiArr","ThiInt"),
                                 coverage_probs = c(0.025, 0.975),
                                 seed = 123) {
        # 共通：DSレベルのモデル妥当性
        P <- build_paths(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i)
        model_chisq <- test_nhpp_chisq(P$counts, P$obs, dt = 1)
        model_rsks  <- test_rescaling_uniform(P$obs, dt = 1)
        
        # 各BSMethodの要約
        rows <- list()
        details <- list()
        
        for (m in methods) {
                paths <- build_paths(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i, BSMethod = m)
                if (!file.exists(paths$ie_each)) {
                        rows[[m]] <- data.frame(
                                method = m,
                                pit_ks_D = NA_real_, pit_ks_p = NA_real_,
                                pit_chisq_p = NA_real_,
                                coverage = NA_real_, k = NA_integer_, n = NA_integer_,
                                cov_binom_p = NA_real_,
                                note = "IEEachPnt file not found",
                                stringsAsFactors = FALSE
                        )
                        next
                }
                # PIT
                pit <- compute_pit_vector(paths$ie_each, paths$obs, seed = seed)
                ks  <- ks_pit_with_jitter(pit)
                chi <- chisq_uniformity(pit, nbins = 10)
                
                # Coverage
                cov <- test_coverage_from_ieeach(paths$ie_each, paths$obs, coverage_probs)
                
                rows[[m]] <- data.frame(
                        method = m,
                        pit_ks_D   = unname(ks$statistic),
                        pit_ks_p   = ks$p.value,
                        pit_chisq_p= chi$p.value,
                        coverage   = cov$coverage,
                        k          = cov$k,
                        n          = cov$n,
                        cov_binom_p= cov$binom_p,
                        note       = "",
                        stringsAsFactors = FALSE
                )
                details[[m]] <- list(pit = pit, ks = ks, chisq = chi, coverage = cov, paths = paths)
        }
        
        out_df <- do.call(rbind, rows)
        rownames(out_df) <- NULL
        list(
                summary = out_df[order(out_df$pit_ks_p), ],
                model_level = list(chisq = model_chisq, rescaling_ks = model_rsks),
                details = details
        )
}
 # ========= 可視化＆表の出力：base R だけで動く =========
 # ---------- タグ生成 ----------
 make_tag <- function(DS, thresholdMode, thresholdName, i, FileName=NULL, dt_code=NULL) {
         parts <- c(
                 sprintf("DS%s", DS),
                 paste0(thresholdMode),
                 paste0(thresholdName),
                 sprintf("i%s", i),
                 if (!is.null(FileName)) paste0("F", FileName) else NULL,
                 if (!is.null(dt_code))  paste0("k", dt_code)  else NULL
         )
         paste(parts[!is.na(parts)], collapse = "_")
 }
 
 # ---------- 置換：保存レンダラ（タグ＆フォルダ対応） ----------
 render_boot_comparison_tagged <- function(res,
                                           out_root = "IE-output/_summaries",
                                           DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i) {
         tag <- make_tag(DS, thresholdMode, thresholdName, i, FileName, dt_code)
         out_dir <- file.path(out_root, GroupDir, FileName, tag)
         if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
         
         df <- res$summary
         rownames(df) <- NULL
         # パラメータ列を付与して保存
         tbl <- data.frame(
                 DS = DS, GroupDir = GroupDir, FileName = FileName,
                 dt_code = dt_code, thresholdMode = thresholdMode, thresholdName = thresholdName, i = i,
                 df, stringsAsFactors = FALSE
         )
         write.csv(tbl, file.path(out_dir, paste0("boot_methods_summary_", tag, ".csv")), row.names = FALSE)
         
         # 図：base R
         # 図：base R（安全ガード付き）
         # 内部ヘルパ：全てNAならプレースホルダPNGを出す
         .plot_bar_safe <- function(file, heights, names, main, ylab,
                                    ylim = NULL, hline = NULL, digits = 3) {
           finite <- is.finite(heights)
           if (!any(finite)) {
             png(file, width = 900, height = 600, res = 120)
             par(mar = c(4, 4, 2, 1))
             plot.new(); title(main = paste0(main, " (no finite data)"))
             dev.off()
             return(invisible(FALSE))
           }
           h <- heights
           # デフォルトylim：0〜(有限値の1.05倍)
           if (is.null(ylim)) {
             ymax <- max(h[finite], na.rm = TRUE)
             ylim <- c(0, if (is.finite(ymax) && ymax > 0) ymax * 1.05 else 1)
           }
           png(file, width = 900, height = 600, res = 120)
           par(mar = c(6, 4, 2, 1))
           bp <- barplot(height = h, names.arg = names, las = 2,
                         ylim = ylim, main = main, ylab = ylab, xlab = "")
           if (!is.null(hline)) abline(h = hline, lty = 2)
           # ラベル：有限値のみ数値、NAは "NA"
           lbl <- ifelse(is.finite(h), round(h, digits), "NA")
           # 棒の上に少し余白を取って表示
           text(x = bp, y = pmin(pmax(h, ylim[1]), ylim[2]), labels = lbl, pos = 3, cex = 0.8)
           dev.off()
           invisible(TRUE)
         }
         
         # --- Coverage（0〜1に固定、名目95%のガイドライン） ---
         .plot_bar_safe(
           file  = file.path(out_dir, paste0("coverage_bar_", tag, ".png")),
           heights = df$coverage,
           names   = df$method,
           main    = sprintf("Coverage by method (%s)", tag),
           ylab    = "Coverage",
           ylim    = c(0, 1),
           hline   = 0.95
         )
         
         # --- PIT-KS -log10 p ---
         ks_logp <- -log10(pmax(df$pit_ks_p, .Machine$double.xmin))
         .plot_bar_safe(
           file  = file.path(out_dir, paste0("pit_ks_logp_bar_", tag, ".png")),
           heights = ks_logp,
           names   = df$method,
           main    = sprintf("PIT-KS -log10 p (%s)", tag),
           ylab    = "-log10(p)",
           hline   = -log10(0.05)
         )
         
         # --- Coverage binom -log10 p ---
         cov_logp <- -log10(pmax(df$cov_binom_p, .Machine$double.xmin))
         .plot_bar_safe(
           file  = file.path(out_dir, paste0("coverage_binom_logp_bar_", tag, ".png")),
           heights = cov_logp,
           names   = df$method,
           main    = sprintf("Coverage binom -log10 p (%s)", tag),
           ylab    = "-log10(p)",
           hline   = -log10(0.05)
         )
         
         sink(file.path(out_dir, paste0("model_level_summary_", tag, ".txt")))
         cat("== Model-level diagnostics ==\n\n")
         print(res$model_level$chisq)
         cat("\n-- Rescaling KS --\n")
         print(res$model_level$rescaling_ks)
         sink(NULL)
         
         invisible(list(
                 out_dir = out_dir, tag = tag,
                 files = list(
                         csv = file.path(out_dir, paste0("boot_methods_summary_", tag, ".csv")),
                         fig_coverage = file.path(out_dir, paste0("coverage_bar_", tag, ".png")),
                         fig_pitks = file.path(out_dir, paste0("pit_ks_logp_bar_", tag, ".png")),
                         fig_covp = file.path(out_dir, paste0("coverage_binom_logp_bar_", tag, ".png")),
                         model_txt = file.path(out_dir, paste0("model_level_summary_", tag, ".txt"))
                 )
         ))
 }
 
 run_and_save_compare <- function(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i,
                                  methods = c("COS","GS-","GS+","OrdSta","ThiArr","ThiInt"),
                                  coverage_probs = c(0.025, 0.975),
                                  seed = 123,
                                  out_root = "IE-output/_summaries") {
         set.seed(seed)
         res <- compare_boot_methods(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i,
                                     methods = methods, coverage_probs = coverage_probs, seed = seed)
         out <- render_boot_comparison_tagged(res, out_root,
                                              DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i)
         message("Saved outputs under: ", out$out_dir)
         res  # 結果オブジェクトも返す（コンソールでも確認可）
 }
# ========== 複数ケースまとめ ==========

collect_compare_results <- function(cases,
                                    methods = c("COS","GS-","GS+","OrdSta","ThiArr","ThiInt"),
                                    coverage_probs = c(0.025,0.975),
                                    seed=123,
                                    out_root="IE-output/_summaries") {
        all_rows <- list()
        for (cc in seq_len(nrow(cases))) {
                DS <- cases$DS[cc]; GroupDir <- cases$GroupDir[cc]; FileName <- cases$FileName[cc]
                dt_code <- cases$dt_code[cc]; thresholdMode <- cases$thresholdMode[cc]
                thresholdName <- cases$thresholdName[cc]; i <- cases$i[cc]
                
                res <- compare_boot_methods(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i,
                                            methods=methods, coverage_probs=coverage_probs, seed=seed)
                
                df <- res$summary
                df$DS <- DS; df$FileName <- FileName; df$dt_code <- dt_code
                df$thresholdMode <- thresholdMode; df$thresholdName <- thresholdName; df$i <- i
                all_rows[[cc]] <- df
        }
        big_df <- do.call(rbind, all_rows)
        # 保存
        out_csv <- file.path(out_root, "all_cases_summary.csv")
        write.csv(big_df, out_csv, row.names=FALSE)
        message("Saved combined summary to: ", out_csv)
        return(big_df)
}

# =========================================================
# まとめ実行 → ケース横断の総括表＆ヒートマップを自動出力
# =========================================================
 save_heatmaps_gg <- function(df, out_dir, tag, fn) {
         .ensure_dir(out_dir)
         df$i_num  <- suppressWarnings(as.numeric(as.character(df$i)))
         df$method <- factor(df$method, levels = .boot_method_order)
         
         # 共通テーマ（白背景）
         white_bg <- theme_minimal(base_size = 14) +
                 theme(
                         panel.background = element_rect(fill="white", colour=NA),
                         plot.background  = element_rect(fill="white", colour=NA)
                 )
         
         # Coverage
         p1 <- ggplot(df, aes(x = method, y = i_num, fill = coverage)) +
                 geom_tile() +
                 scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
                 scale_y_continuous(breaks = sort(unique(df$i_num))) +
                 labs(title = sprintf("Coverage (%s_%s)", tag, fn), x = "method", y = "i") +
                 white_bg
         ggsave(file.path(out_dir, sprintf("coverage_heatmap_gg_%s_%s.png", tag, fn)),
                p1, width = 6, height = 5, dpi = 120, bg="white")
         
         # PIT-KS -log10 p
         df$pit_ks_logp <- -log10(pmax(df$pit_ks_p, .Machine$double.xmin))
         p2 <- ggplot(df, aes(x = method, y = i_num, fill = pit_ks_logp)) +
                 geom_tile() +
                 scale_fill_gradient(low = "white", high = "red") +
                 scale_y_continuous(breaks = sort(unique(df$i_num))) +
                 labs(title = sprintf("PIT-KS -log10 p (%s_%s)", tag, fn), x = "method", y = "i") +
                 white_bg
         ggsave(file.path(out_dir, sprintf("pitks_heatmap_gg_%s_%s.png", tag, fn)),
                p2, width = 6, height = 5, dpi = 120, bg="white")
         
         # Coverage binom -log10 p
         df$covbinom_logp <- -log10(pmax(df$cov_binom_p, .Machine$double.xmin))
         p3 <- ggplot(df, aes(x = method, y = i_num, fill = covbinom_logp)) +
                 geom_tile() +
                 scale_fill_gradient(low = "white", high = "darkgreen") +
                 scale_y_continuous(breaks = sort(unique(df$i_num))) +
                 labs(title = sprintf("Coverage binom -log10 p (%s_%s)", tag, fn), x = "method", y = "i") +
                 white_bg
         ggsave(file.path(out_dir, sprintf("covbinom_heatmap_gg_%s_%s.png", tag, fn)),
                p3, width = 6, height = 5, dpi = 120, bg="white")
         
         # AD
         df$nhpp_ad_logp <- -log10(pmax(df$nhpp_ad_p, .Machine$double.xmin))
         p4 <- ggplot(df, aes(x=method, y=i_num, fill=nhpp_ad_logp)) +
                 geom_tile() +
                 scale_fill_gradient(low="white", high="purple") +
                 labs(title=sprintf("NHPP AD -log10 p (%s_%s)", tag, fn), x="method", y="i") +
                 white_bg
         ggsave(file.path(out_dir, sprintf("nhpp_ad_heatmap_gg_%s_%s.png", tag, fn)),
                p4, width=6, height=5, dpi=120, bg="white")
         
         # CvM
         df$nhpp_cvm_logp <- -log10(pmax(df$nhpp_cvm_p, .Machine$double.xmin))
         p5 <- ggplot(df, aes(x=method, y=i_num, fill=nhpp_cvm_logp)) +
                 geom_tile() +
                 scale_fill_gradient(low="white", high="brown") +
                 labs(title=sprintf("NHPP CvM -log10 p (%s_%s)", tag, fn), x="method", y="i") +
                 white_bg
         ggsave(file.path(out_dir, sprintf("nhpp_cvm_heatmap_gg_%s_%s.png", tag, fn)),
                p5, width=6, height=5, dpi=120, bg="white")
 }

 # 返り値: p値のみを返す（表作成しやすい）
 nhpp_diagnostics <- function(path_counts, path_lambda, dt = 1, min_exp = 5, jitter_eps = 1e-12) {
         read_counts_vec <- function(path_counts) {
                 df <- utils::read.table(path_counts, header = FALSE, sep = "",
                                         check.names = FALSE, comment.char = "",
                                         quote = "", stringsAsFactors = FALSE)
                 as.numeric(df[[2]])
         }
         read_lambda_vec <- function(path_lambda) {
                 df <- utils::read.table(path_lambda, header = TRUE, sep = ",",
                                         check.names = FALSE, comment.char = "",
                                         stringsAsFactors = FALSE)
                 as.numeric(df[[1]])
         }
         
         counts <- read_counts_vec(path_counts)
         lam    <- read_lambda_vec(path_lambda)
         stopifnot(length(counts) == length(lam))
         
         # χ²（参考: p値だけ返す）
         exp_counts <- lam * dt
         O <- c(); E <- c(); accO <- 0; accE <- 0
         for (j in seq_along(counts)) {
                 accO <- accO + counts[j]; accE <- accE + exp_counts[j]
                 if (accE >= min_exp || j == length(counts)) {
                         O <- c(O, accO); E <- c(E, accE); accO <- 0; accE <- 0
                 }
         }
         df_chi <- length(O) - 1
         chisq  <- sum((O - E)^2 / pmax(E, .Machine$double.eps))
         p_chi  <- stats::pchisq(chisq, df = df_chi, lower.tail = FALSE)
         
         # rescaling → U(0,1) になるはず
         u <- 1 - exp(-(lam * dt))
         u <- pmin(pmax(u + stats::runif(length(u), -jitter_eps, jitter_eps), 0), 1)
         u <- u[is.finite(u) & u > 0 & u < 1]
         
         ks_res  <- stats::ks.test(u, "punif")
         ad_res  <- tryCatch(goftest::ad.test(u, null = "punif"),  error = function(e) NULL)
         cvm_res <- tryCatch(goftest::cvm.test(u, null = "punif"), error = function(e) NULL)
         
         list(
                 nhpp_chisq_p = p_chi,
                 nhpp_ks_p    = ks_res$p.value,
                 nhpp_ad_p    = if (!is.null(ad_res))  ad_res$p.value  else NA_real_,
                 nhpp_cvm_p   = if (!is.null(cvm_res)) cvm_res$p.value else NA_real_
         )
 }
 
 # ================= NHPP ヒートマップ（KS/AD/CvM の −log10 p） =================
 save_nhpp_heatmap_gg <- function(nhpp_df_unique_i, out_dir, tag, fn) {
         # nhpp_df_unique_i: i ごとに 1 行（nhpp_ks_p, nhpp_ad_p, nhpp_cvm_p を列にもつ）
         .ensure_dir(out_dir)
         if (is.null(nhpp_df_unique_i) || !nrow(nhpp_df_unique_i)) {
                 message("save_nhpp_heatmap_gg: empty input; skip"); return(invisible(FALSE))
         }
         
         # 必要列チェック
         need <- c("i","nhpp_ks_p","nhpp_ad_p","nhpp_cvm_p")
         miss <- setdiff(need, names(nhpp_df_unique_i))
         if (length(miss)) stop("save_nhpp_heatmap_gg: missing columns: ", paste(miss, collapse=", "))
         
         df_long <- reshape(
                 data = nhpp_df_unique_i[, need],
                 varying = list(2:4),
                 v.names = "pval",
                 timevar = "test",
                 times = c("KS", "AD", "CvM"),
                 direction = "long"
         )
         df_long$logp  <- -log10(pmax(df_long$pval, .Machine$double.xmin))
         df_long$i_num <- suppressWarnings(as.numeric(as.character(df_long$i)))
         df_long$test  <- factor(df_long$test, levels = c("KS","AD","CvM"))
         
         # i を数値化できなかった行は除外
         df_long <- df_long[is.finite(df_long$i_num), , drop = FALSE]
         if (!nrow(df_long)) {
                 message("save_nhpp_heatmap_gg: no finite rows after parsing; skip"); return(invisible(FALSE))
         }
         
         # 白背景テーマ
         white_bg <- ggplot2::theme_minimal(base_size = 14) +
                 ggplot2::theme(
                         panel.background = ggplot2::element_rect(fill = "white", colour = NA),
                         plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
                 )
         
         p <- ggplot2::ggplot(df_long, ggplot2::aes(x = test, y = i_num, fill = logp)) +
                 ggplot2::geom_tile() +
                 ggplot2::scale_fill_gradient(low = "white", high = "purple") +
                 ggplot2::scale_y_continuous(breaks = sort(unique(df_long$i_num))) +
                 ggplot2::labs(
                         title = sprintf("NHPP diagnostics -log10(p) (%s_%s)", tag, fn),
                         x = "test", y = "i", fill = "-log10(p)"
                 ) +
                 white_bg
         
         ggplot2::ggsave(
                 file.path(out_dir, sprintf("nhpp_diag_heatmap_gg_%s_%s.png", tag, fn)),
                 p, width = 6, height = 5, dpi = 120, bg = "white"
         )
         invisible(TRUE)
 }
 
 
# 1) すべてのケースを回して1つの大テーブルにまとめる（ついでにケース別の個別出力も実施）
 run_many_and_save <- function(cases,
                               methods = .boot_method_order,
                               coverage_probs = c(0.025, 0.975),
                               seed = 123,
                               out_root = "IE-output/_summaries") {
         set.seed(seed)
         big_rows <- list()
         
         for (cc in seq_len(nrow(cases))) {
                 DS <- as.character(cases$DS[cc])
                 GroupDir <- as.character(cases$GroupDir[cc])
                 FileName <- as.character(cases$FileName[cc])
                 dt_code <- as.character(cases$dt_code[cc])
                 thresholdMode <- as.character(cases$thresholdMode[cc])
                 thresholdName <- as.character(cases$thresholdName[cc])
                 i <- as.character(cases$i[cc])
                 
                 # 1) 各法の比較（coverage / PIT など）
                 res <- compare_boot_methods(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i,
                                             methods = methods, coverage_probs = coverage_probs, seed = seed)
                 
                 # 2) NHPP 診断（KS/AD/CvM）— このケースの「観測」ファイルを直接使う
                 P <- build_paths(DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i)
                 PATH_COUNTS <- P$counts
                 PATH_OBS    <- P$obs
                 diag <- nhpp_diagnostics(PATH_COUNTS, PATH_OBS, dt = 1)
                 
                 # 3) ケース個別の保存（棒グラフはこの関数の中で出力される想定）
                 render_boot_comparison_tagged(res, out_root,
                                               DS, GroupDir, FileName, dt_code, thresholdMode, thresholdName, i)
                 
                 # 4) 総括用に集約（各法の行に NHPP の p値を同じ値で付与）
                 df <- res$summary
                 df$DS <- DS; df$GroupDir <- GroupDir; df$FileName <- FileName
                 df$dt_code <- dt_code; df$thresholdMode <- thresholdMode; df$thresholdName <- thresholdName; df$i <- i
                 df$nhpp_chisq_p <- diag$nhpp_chisq_p
                 df$nhpp_ks_p    <- diag$nhpp_ks_p
                 df$nhpp_ad_p    <- diag$nhpp_ad_p
                 df$nhpp_cvm_p   <- diag$nhpp_cvm_p
                 big_rows[[length(big_rows) + 1]] <- df
         }
         
         big_df <- do.call(rbind, big_rows)
         
         # === 出力先は <GroupDir>/<FileName> 配下。CSVは上書き保存 ===
         combos <- unique(big_df[, c("GroupDir", "FileName")])
         for (r in seq_len(nrow(combos))) {
                 gd <- combos$GroupDir[r]; fn <- combos$FileName[r]
                 out_dir <- file.path(out_root, gd, fn)
                 .ensure_dir(out_dir)
                 
                 part <- subset(big_df, GroupDir == gd & FileName == fn)
                 write.csv(part, file.path(out_dir, "all_cases_summary.csv"), row.names = FALSE)
                 
                 # DS × 閾値 × dt_code ごとにヒートマップを保存（ggplot2）
                 combos2 <- unique(part[, c("DS","thresholdMode","thresholdName","dt_code")])
                 for (s in seq_len(nrow(combos2))) {
                         DSv <- combos2$DS[s]; tm <- combos2$thresholdMode[s]
                         tn  <- combos2$thresholdName[s]; dk <- combos2$dt_code[s]
                         tag <- sprintf("DS%s_%s_%s_k%s", DSv, tm, tn, dk)
                         
                         sub <- subset(part, DS == DSv & thresholdMode == tm & thresholdName == tn & dt_code == dk)
                         sub$method <- factor(sub$method, levels = .boot_method_order)
                         sub$i_num  <- suppressWarnings(as.numeric(as.character(sub$i)))
                         
                         # 既存3種（Coverage / PIT-KS / Coverage-binomial）
                         save_heatmaps_gg(sub, out_dir, tag, fn)
                         
                         # NHPP（KS/AD/CvM）— iごとに1行へ縮約しつつヒートマップ保存
                         nhpp_i <- unique(sub[, c("i", "nhpp_ks_p", "nhpp_ad_p", "nhpp_cvm_p")])
                         save_nhpp_heatmap_gg(nhpp_i, out_dir, tag, fn)
                 }
                 message("Saved combined summary & ggplot2 heatmaps under: ", out_dir)
         }
         
         invisible(big_df)
 }
 

#############################
# 使い方（例）
#############################
# 例）あなたのケース
# res <- compare_boot_methods(
#   DS = "1", GroupDir = "DT_Ans_WSE", FileName = "A1",
#   dt_code = "A1", thresholdMode = "h", thresholdName = "ut", i = 2
# )

# 結果の見方
# res$summary        # ← 各法のPIT-KS, PIT-χ², Coverage, Binom p が並ぶ要約
# res$model_level    # ← DSレベルの χ² と rescaling KS（全法共通）
# res$details$COS$ks # ← 任意法の詳細（PITベクトルなど）

# all_res <- collect_compare_results(cases)
# ggplot(all_res, aes(x=method, y=i, fill=coverage)) +
#         geom_tile() + scale_fill_gradient2(mid=0.95, limits=c(0,1)) +
#         facet_wrap(~DS)

 for (ii in c(2,3,4,5,6)) {
         run_and_save_compare(DS="4", GroupDir="NDT_WSE", FileName="",
                              dt_code="H", thresholdMode="h", thresholdName="ldt", i=ii)
 }
 
 # DS = 1/4, i = 2:5/6 を例に（必要に応じて調整）
 cases <- expand.grid(
         #DS = as.character(1:4),
         DS="4",
         GroupDir = "NDT_WSE",
         FileName = "",         # 例：A1（他があればここをベクトル化） ""だとHやFiにできる
         dt_code = "H",
         thresholdMode = "h",
         thresholdName = "ldt",
         i = 2:6,
         stringsAsFactors = FALSE
 )

all_df <- run_many_and_save(cases)

#出力結果についての考察方法は、詳しくはkentei考察.txt