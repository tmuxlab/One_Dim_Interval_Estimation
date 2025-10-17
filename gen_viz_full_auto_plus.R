library(ggplot2)
library(reshape2)

source("run_PEIE_addons_generic.R")  # 既存（生成器）
source("run_PEIE_full_auto_plus.R")  # いま貼った拡張版

plot_peie_heatmap <- function(summary_csv, metric=c("KS_PercAns_p","KS_tQuick_p","Cov_PercAns_p","Cov_tQuick_p")){
        if (!file.exists(summary_csv)) stop("Not found: ", summary_csv)
        df <- read.csv(summary_csv, check.names=FALSE, stringsAsFactors=FALSE)
        
        for (m in metric){
                if (!m %in% names(df)) {
                        warning("Metric ", m, " not in file, skip."); next
                }
                dat <- df[,c("COS","J",m)]
                colnames(dat) <- c("COS","J","pval")
                dat$flag <- cut(dat$pval,
                                breaks=c(-Inf,0.05,0.1,1),
                                labels=c("p<0.05","0.05≤p<0.1","p≥0.1"))
                dat$J <- factor(dat$J, levels=sort(unique(dat$J)))
                p <- ggplot(dat, aes(x=J, y=COS, fill=flag))+
                        geom_tile(color="white")+
                        geom_text(aes(label=sprintf("%.3f", pval)), size=3)+
                        scale_fill_manual(values=c("p<0.05"="red","0.05≤p<0.1"="orange","p≥0.1"="green"))+
                        labs(title=paste("Heatmap:", m), x="J", y="COS")+
                        theme_minimal()+
                        theme(axis.text.x=element_text(angle=0,hjust=0.5))
                print(p)
                
                # 保存
                outfile <- sub("\\.csv$", paste0("-HEATMAP-",m,".png"), summary_csv)
                ggsave(outfile, plot=p, width=6, height=4, dpi=150)
                message("Saved heatmap: ", outfile)
        }
}

#PIT系×Coverage系のマルチパネル図
plot_peie_multipanel <- function(summary_csv, png_out = NULL){
        stopifnot(file.exists(summary_csv))
        suppressWarnings(suppressMessages(require(ggplot2)))
        suppressWarnings(suppressMessages(require(gridExtra)))
        
        df <- read.csv(summary_csv, check.names = FALSE, stringsAsFactors = FALSE)
        
        build_heat <- function(metric, title){
                if (!metric %in% names(df)) return(ggplot() + theme_void() + ggtitle(paste(title, "(not found)")))
                dat <- df[, c("COS","J", metric)]
                colnames(dat) <- c("COS","J","pval")
                dat$flag <- cut(dat$pval, breaks = c(-Inf, 0.05, 0.10, 1),
                                labels = c("p<0.05", "0.05≤p<0.1", "p≥0.1"))
                dat$J   <- factor(dat$J, levels = sort(unique(dat$J)))
                dat$COS <- factor(dat$COS, levels = rev(sort(unique(dat$COS))))  # 上からC0S…順
                ggplot(dat, aes(x = J, y = COS, fill = flag)) +
                        geom_tile(color = "white", linewidth = .4) +
                        geom_text(aes(label = ifelse(is.na(pval), "", sprintf("%.3f", pval))), size = 3) +
                        scale_fill_manual(values = c("p<0.05" = "red", "0.05≤p<0.1" = "orange", "p≥0.1" = "green"),
                                          na.value = "grey80") +
                        labs(title = title, x = "J", y = "COS") +
                        theme_minimal(base_size = 12) +
                        theme(legend.position = "right",
                              panel.grid = element_blank())
        }
        
        g1 <- build_heat("KS_PercAns_p",  "PIT (KS) — PercAns")
        g2 <- build_heat("KS_tQuick_p",   "PIT (KS) — tQuick")
        g3 <- build_heat("Cov_PercAns_p", "Coverage (Binom) — PercAns")
        g4 <- build_heat("Cov_tQuick_p",  "Coverage (Binom) — tQuick")
        
        if (is.null(png_out)) {
                png_out <- sub("\\.csv$", "-HEATMAP-4in1.png", summary_csv)
        }
        g <- grid.arrange(g1, g2, g3, g4, ncol = 2)
        ggsave(png_out, g, width = 12, height = 8, dpi = 150)
        message("Saved: ", png_out)
        invisible(png_out)
}

compare_ci_width_position <- function(
                DS, thrMode, thrName, J,
                method_dir, variant = NULL,     # NDT_WSE は variant=NULL、DT_* は A1/A2/A3 or B1/B2
                cos_tags = c("COS","GS-","GS+","OrdSta","ThiArr","ThiInt")
){
        # ---- パス&タグ ----
        is_ndt <- identical(method_dir, "NDT_WSE")
        if (!is_ndt && thrName == "ldt") stop("thrName='ldt' は NDT_WSE 専用です。")
        kcode_eff <- if (is_ndt) "kH" else paste0("k", variant)
        base_dir  <- if (is_ndt) file.path("IE-output", method_dir) else file.path("IE-output", method_dir, variant)
        
        .path <- function(tag, cos_tag){
                file.path(base_dir, sprintf("DS%s-IEEachPnt-%s-%s-%s-%s-%s(J=%s).csv",
                                            DS, tag, cos_tag, kcode_eff, thrMode, thrName, J))
        }
        
        # ---- CI読込ヘルパ ----
        read_ci <- function(tag, cos_tag){
                f <- .path(tag, cos_tag)
                if (!file.exists(f)) return(NULL)
                df <- read.csv(f, check.names = FALSE)
                if (!all(c("lo","hi") %in% names(df))) return(NULL)
                df
        }
        
        # ---- 要約1本（あるCOSに対して） ----
        summarize_one <- function(cos_tag){
                A <- read_ci("PercAns", cos_tag)
                TQ <- read_ci("tQuick", cos_tag)
                if (is.null(A) || is.null(TQ)) return(NULL)
                
                widthA <- pmax(0, A$hi - A$lo)
                widthT <- pmax(0, TQ$hi - TQ$lo)
                midA   <- (A$hi + A$lo)/2
                midT   <- (TQ$hi + TQ$lo)/2
                
                # 重なり率（overlap/union）
                overlap <- pmax(0, pmin(A$hi, TQ$hi) - pmax(A$lo, TQ$lo))
                union   <- pmax(widthA + widthT - overlap, 1e-12)
                jacc    <- overlap / union
                
                data.frame(
                        DS = DS, Method = method_dir, Variant = if (is_ndt) "H" else variant,
                        COS = cos_tag, J = J, thrMode = thrMode, thrName = thrName,
                        # 幅の指標
                        mean_width_PercAns = mean(widthA, na.rm = TRUE),
                        mean_width_tQuick  = mean(widthT, na.rm = TRUE),
                        mean_width_diff    = mean(widthT - widthA, na.rm = TRUE),
                        prop_tQuick_wider  = mean(widthT > widthA, na.rm = TRUE),
                        rel_widen_mean     = mean((widthT / pmax(widthA, 1e-12)) - 1, na.rm = TRUE),
                        # 位置の指標
                        mean_mid_absdiff   = mean(abs(midT - midA), na.rm = TRUE),
                        rmse_mid_diff      = sqrt(mean((midT - midA)^2, na.rm = TRUE)),
                        # 区間の重なり
                        mean_overlap_ratio = mean(jacc, na.rm = TRUE)
                )
        }
        
        rows <- lapply(cos_tags, summarize_one)
        out  <- do.call(rbind, rows)
        
        # 保存
        out_file <- file.path(base_dir, sprintf("DS%s-CI-COMPARE-%s-%s-%s(J=%s).csv", DS, kcode_eff, thrMode, thrName, J))
        write.csv(out, out_file, row.names = FALSE)
        message("Saved CI comparison: ", out_file)
        out
}

# 例：DS1、thrMode="h"、J=2:5、体制 = NDT_H + A3 + B2、COSタグ=既定6種
run_PEIE_all(
        DS=4, thrMode="h",
        J_vec=2:6,
        tasks=list(
                list(method_dir="NDT_WSE",    variant=NULL, thrName="ldt"),
                list(method_dir="DT_Ans_WSE", variant="A3", thrName="ut"),
                list(method_dir="DT_Bar_WSE", variant="B2", thrName="ut")
        )
)

# NDT_WSE のサマリーを可視化
#plot_peie_heatmap("IE-output/NDT_WSE/DS1-GRID-SUMMARY-NDT_WSE(thrMode=h).csv")
plot_peie_multipanel("IE-output/NDT_WSE/DS4-GRID-SUMMARY-NDT_WSE(thrMode=h).csv")
plot_peie_multipanel("IE-output/DT_Ans_WSE/A3/DS4-GRID-SUMMARY-DT_Ans_WSE(thrMode=h).csv")
plot_peie_multipanel("IE-output/DT_Bar_WSE/B2/DS4-GRID-SUMMARY-DT_Bar_WSE(thrMode=h).csv")

#ここからはPercAns と tQuickの比較、同じでないことがわかる。
# # 例: NDT_WSE / DS1 / thrMode="h" / thrName="ldt" / J=2
# cmp <- compare_ci_width_position(
#         DS=1, thrMode="h", thrName="ldt", J=2,
#         method_dir="NDT_WSE", variant=NULL
# )
# print(cmp)
# 
# # 例: DT_Ans_WSE / A3
# compare_ci_width_position(
#         DS=1, thrMode="h", thrName="ut", J=2,
#         method_dir="DT_Ans_WSE", variant="A3"
# )


# # ================== ハイライトスクリプト ==================
# highlight_peie_summary <- function(
#                 summary_csv,
#                 p_alpha = 0.05,   # 有意判定（強）
#                 p_warn  = 0.10    # 注意喚起（弱）
# ){
#         if (!file.exists(summary_csv)) stop("Not found: ", summary_csv)
#         x <- read.csv(summary_csv, check.names = FALSE, stringsAsFactors = FALSE)
#         
#         # 補助: 安全に NA を FALSE 扱いするヘルパ
#         lt <- function(a, thr) ifelse(is.na(a), FALSE, a < thr)
#         
#         # 主要フラグ
#         x$PIT_bad_strict  <- lt(x$KS_PercAns_p, p_alpha) | lt(x$KS_tQuick_p, p_alpha)
#         x$PIT_bad_warn    <- lt(x$KS_PercAns_p, p_warn)  | lt(x$KS_tQuick_p, p_warn)
#         
#         x$COV_bad_strict  <- lt(x$Cov_PercAns_p, p_alpha) | lt(x$Cov_tQuick_p, p_alpha)
#         x$COV_bad_warn    <- lt(x$Cov_PercAns_p, p_warn)  | lt(x$Cov_tQuick_p, p_warn)
#         
#         # 併せて、どちらで引っかかったかを文字で持つ（最終CSVに入れる）
#         flag_reason <- function(row){
#                 rs <- c()
#                 if (!is.na(row["KS_PercAns_p"]) && row["KS_PercAns_p"] < p_alpha) rs <- c(rs, "PIT(PercAns)")
#                 if (!is.na(row["KS_tQuick_p"])  && row["KS_tQuick_p"]  < p_alpha) rs <- c(rs, "PIT(tQuick)")
#                 if (!is.na(row["Cov_PercAns_p"])&& row["Cov_PercAns_p"]< p_alpha) rs <- c(rs, "Cov(PercAns)")
#                 if (!is.na(row["Cov_tQuick_p"]) && row["Cov_tQuick_p"] < p_alpha) rs <- c(rs, "Cov(tQuick)")
#                 if (length(rs)==0) "" else paste(rs, collapse = "; ")
#         }
#         x$Strict_Flags <- apply(x, 1, flag_reason)
#         
#         # コンソール向けの簡易ハイライト表（COS×J×体制）
#         # 記号: "!!"=strict, "!"=warn, ""=OK
#         x$PIT_mark <- ifelse(x$PIT_bad_strict, "!!", ifelse(x$PIT_bad_warn, "!", ""))
#         x$COV_mark <- ifelse(x$COV_bad_strict, "!!", ifelse(x$COV_bad_warn, "!", ""))
#         
#         show_cols <- c("Method","Variant","COS","J","thrMode","thrName",
#                        "KS_PercAns_p","KS_tQuick_p","PIT_mark",
#                        "Cov_PercAns_prop","Cov_tQuick_prop","Cov_PercAns_p","Cov_tQuick_p","COV_mark")
#         show_cols <- show_cols[show_cols %in% names(x)]
#         x_show <- x[, show_cols, drop = FALSE]
#         
#         # 並べ替え（Method→Variant→COS→J）
#         ord <- order(x$Method, x$Variant, x$COS, x$J)
#         x_show <- x_show[ord, ]
#         
#         # 出力（コンソール）
#         cat("\n=== HIGHLIGHTS (", basename(summary_csv), ") ===\n", sep = "")
#         print(x_show, row.names = FALSE)
#         cat("\nLegend: PIT_mark/COV_mark — '!!' strict (p<", p_alpha,
#             "), '!' warn (p<", p_warn, "), '' ok\n\n", sep = "")
#         
#         # 重要列だけ CSV として保存
#         out_csv <- sub("\\.csv$", "-HIGHLIGHTS.csv", summary_csv)
#         write.csv(cbind(x_show, Strict_Flags = x$Strict_Flags[ord]), out_csv, row.names = FALSE)
#         message("Saved highlights: ", out_csv)
#         
#         invisible(list(table = x_show, path = out_csv))
# }
# # ================== 使い方例 ==================
# # 例：NDT_WSE のグリッドサマリー
# highlight_peie_summary("IE-output/NDT_WSE/DS1-GRID-SUMMARY-NDT_WSE(thrMode=h).csv")
