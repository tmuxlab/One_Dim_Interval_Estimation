collect_bca_coverage_grid <- function(
                DS, method_dir, variant = NULL,
                thrMode, thrName,
                J_vec,
                cos_tags = c("COS","GS-","GS+","OrdSta","ThiArr","ThiInt"),
                make_wide = TRUE,       # TRUE なら prop を Wide でも保存
                verbose = TRUE,
                stop_if_missing = FALSE # TRUE なら欠損時に stop
){
        # 体制ごとのパス/kcode
        is_ndt <- identical(method_dir, "NDT_WSE")
        if (!is_ndt && identical(thrName, "ldt")) {
                stop("thrName='ldt' は NDT_WSE 専用です（DT系では不可）。")
        }
        kcode <- if (is_ndt) "kH" else paste0("k", variant)
        base_ie <- if (is_ndt) file.path("IE-output", method_dir)
        else        file.path("IE-output", method_dir, variant)
        dir.create(base_ie, recursive = TRUE, showWarnings = FALSE)
        
        # 収集
        rows <- list(); miss <- character(0)
        for (ct in cos_tags) {
                for (JJ in J_vec) {
                        path <- file.path(
                                base_ie,
                                sprintf("DS%s-BCACOMP-%s-%s-%s-%s(J=%s).csv",
                                        DS, ct, kcode, thrMode, thrName, JJ)
                        )
                        if (!file.exists(path)) {
                                msg <- paste0("missing: ", path)
                                miss <- c(miss, msg)
                                if (verbose) message(msg)
                                if (isTRUE(stop_if_missing)) stop(msg)
                                # スキップ
                        } else {
                                df <- try(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
                                if (inherits(df, "try-error") || !all(c("BCA_coverage_inside","BCA_coverage_prop","BCA_coverage_p") %in% names(df))) {
                                        msg <- paste0("invalid file (columns not found): ", path)
                                        miss <- c(miss, msg)
                                        if (verbose) message(msg)
                                        if (isTRUE(stop_if_missing)) stop(msg)
                                } else {
                                        # 1行想定、必要列だけ抽出
                                        take <- c("DS","GroupDir","FileName","COS","J","BSMethod","dt","thrMode","thrName",
                                                  "BCA_coverage_inside","BCA_coverage_prop","BCA_coverage_p")
                                        take <- take[take %in% names(df)]
                                        df2 <- df[, take, drop = FALSE]
                                        rows[[length(rows)+1]] <- df2
                                }
                        }
                }
        }
        
        if (!length(rows)) stop("収集できる BCACOMP がありません。")
        
        # 縦持ちで結合
        out_long <- do.call(rbind, rows)
        
        # 保存（long版）
        out_long_file <- file.path(base_ie,
                                   sprintf("DS%s-BCACOV-GRID-%s(thrMode=%s).csv", DS, method_dir, thrMode)
        )
        write.csv(out_long, out_long_file, row.names = FALSE)
        if (verbose) message("Saved long: ", out_long_file)
        
        # Wide 版（propのみ：COS×J 行列）
        if (isTRUE(make_wide)) {
                # safety: 型を明示
                out_long$COS <- as.character(out_long$COS)
                out_long$J   <- as.integer(out_long$J)
                out_long$BCA_coverage_prop <- as.numeric(out_long$BCA_coverage_prop)
                
                # spread
                J_levels <- sort(unique(out_long$J))
                COS_levels <- unique(out_long$COS)
                # 体裁：COS（行）× J（列）
                wide_prop <- do.call(rbind, lapply(COS_levels, function(ct) {
                        row <- numeric(length(J_levels))
                        for (i in seq_along(J_levels)) {
                                jj <- J_levels[i]
                                v <- out_long$BCA_coverage_prop[out_long$COS==ct & out_long$J==jj]
                                row[i] <- if (length(v)) v[1] else NA_real_
                        }
                        data.frame(COS = ct, setNames(as.list(row), paste0("J", J_levels)), check.names = FALSE)
                }))
                out_wide_file <- file.path(base_ie,
                                           sprintf("DS%s-BCACOV-WIDEprop-%s(thrMode=%s).csv", DS, method_dir, thrMode)
                )
                write.csv(wide_prop, out_wide_file, row.names = FALSE)
                if (verbose) message("Saved wide(prop): ", out_wide_file)
        }
        
        invisible(list(long = out_long_file, wide_prop = if (isTRUE(make_wide)) out_wide_file else NULL,
                       missing = miss))
}

# # 例：DT_Bar_WSE / B2, DS=1, thrMode="h", thrName="ut", J=2:5, 6つのCOS
# collect_bca_coverage_grid(
#         DS=4, method_dir="DT_Bar_WSE", variant="B2",
#         thrMode="h", thrName="ut",
#         J_vec=2:6,
#         cos_tags=c("COS","GS-","GS+","OrdSta","ThiArr","ThiInt"),
#         make_wide=TRUE  # ← ヒートマップ用の wide も同時保存（propのみ）
# )
collect_bca_coverage_grid(
        DS=1, method_dir="DT_Ans_WSE", variant="A3",
        thrMode="h", thrName="ut",
        J_vec=2:5,
        cos_tags=c("COS","GS-","GS+","OrdSta","ThiArr","ThiInt"),
        make_wide=TRUE  # ← ヒートマップ用の wide も同時保存（propのみ）
)
# collect_bca_coverage_grid(
#         DS=4, method_dir="NDT_WSE", variant="H",
#         thrMode="h", thrName="ldt",
#         J_vec=2:6,
#         cos_tags=c("COS","GS-","GS+","OrdSta","ThiArr","ThiInt"),
#         make_wide=TRUE  # ← ヒートマップ用の wide も同時保存（propのみ）
# )
