# ============================================================
# 汎用スタンドアロン：NDT_WSE / DT_Ans_WSE / DT_Bar_WSE 用
# 規約：
#  - NDT_WSE 配下はサブフォルダなし（Hのみ格納）
#  - DT_Ans_WSE / DT_Bar_WSE は variant サブフォルダ（A1/A2/A3, B1/B2）
#  - thrName="ldt" は NDT_WSE 専用
#  - kcode は出力側で "k" + variant（NDTは kH、DTは kA1/kA2/kA3/kB1/kB2）
#  - BSCount <TAG>：NDTは H 固定、DTは variant 固定（フォールバック無し）
# 出力先：IE-output/<method_dir>/<variant?>/DS*-IEEachPnt-...csv, IEBand-...csv
# ============================================================

run_PEIE_addons_generic <- function(
                DS, thrMode, thrName, J,
                method_dir = "NDT_WSE",   # "NDT_WSE" | "DT_Ans_WSE" | "DT_Bar_WSE"
                variant    = NULL,        # NDT: NULL/""/ignored（内部で H 扱い）; DT: "A1"/"A2"/"A3" or "B1"/"B2"
                cos_tag    = "COS",
                alpha = 0.05,
                delta_t = 1,              # 等間隔なら 1。不等間隔なら長さTのベクトル
                transforms = c("Anscombe"),
                emit_tquick = TRUE,
                emit_band_sup = TRUE,
                emit_band_supt = FALSE
){
        # ---------- 体制チェック ----------
        is_ndt <- identical(method_dir, "NDT_WSE")
        is_ans <- identical(method_dir, "DT_Ans_WSE")
        is_bar <- identical(method_dir, "DT_Bar_WSE")
        if (!(is_ndt || is_ans || is_bar)) stop("method_dir は 'NDT_WSE'/'DT_Ans_WSE'/'DT_Bar_WSE' のいずれか。")
        
        # thrName=ldt の制約
        if (!is_ndt && identical(thrName, "ldt")) {
                stop("thrName='ldt' は NDT_WSE 専用です（DT系では使用不可）。")
        }
        
        # variant/kcode/BSCountタグの決定
        if (is_ndt) {
                variant_eff <- "H"               # NDTは常にHだが、サブフォルダは作らない
                kcode_eff   <- "kH"
                bscount_tag <- "H"
        } else {
                if (is.null(variant) || !nzchar(variant)) {
                        stop("DT系では variant を指定してください（DT_Ans_WSE: A1/A2/A3, DT_Bar_WSE: B1/B2）")
                }
                # variantの妥当性
                if (is_ans && !variant %in% c("A1","A2","A3")) stop("DT_Ans_WSE の variant は A1/A2/A3 のみ。")
                if (is_bar && !variant %in% c("B1","B2"))       stop("DT_Bar_WSE の variant は B1/B2 のみ。")
                variant_eff <- variant
                kcode_eff   <- paste0("k", variant_eff)
                bscount_tag <- variant_eff
        }
        
        # ---------- パスユーティリティ ----------
        .p <- function(...) {
                parts <- list(...)
                parts <- parts[!vapply(parts, function(x) is.null(x) || (is.character(x) && x==""), logical(1))]
                do.call(file.path, parts)
        }
        # NDTはサブフォルダなし、DTは variant サブフォルダあり
        .stage_dir <- function(stage){
                if (is_ndt) .p(stage, method_dir) else .p(stage, method_dir, variant_eff)
        }
        
        # ---------- 読み込み ----------
        .read_lambda_hat <- function(){
                f <- .p(.stage_dir("output"),
                        sprintf("DS%s-%s-%s-%s(J=%s).txt", DS, kcode_eff, thrMode, thrName, J))
                if (!file.exists(f)) stop("Not found: ", f)
                df <- try(read.csv(f, header=TRUE, check.names=FALSE, strip.white=TRUE,
                                   stringsAsFactors=FALSE), silent=TRUE)
                if (inherits(df, "try-error")) {
                        df <- try(read.table(f, header=TRUE, check.names=FALSE, strip.white=TRUE,
                                             stringsAsFactors=FALSE), silent=TRUE)
                }
                if (inherits(df, "try-error")) df <- read.table(f, header=FALSE, strip.white=TRUE)
                as.numeric(df[[1]])
        }
        
        .read_BSCount_strict <- function(){
                f <- .p(.stage_dir("BS-output"),
                        sprintf("DS%s-BSCount-%s-%s-%s-%s(J=%s).dat",
                                DS, cos_tag, bscount_tag, thrMode, thrName, J))
                if (!file.exists(f)) stop("Not found: ", f)
                df <- try(read.table(f, header=TRUE, check.names=FALSE, strip.white=TRUE,
                                     stringsAsFactors=FALSE), silent=TRUE)
                if (inherits(df, "try-error")) {
                        df <- read.table(f, header=FALSE, check.names=FALSE, strip.white=TRUE)
                }
                df[] <- lapply(df, function(col) suppressWarnings(as.numeric(col)))
                as.matrix(df)
        }
        
        .ensure_TxB <- function(mat, T_expected){
                # B×T 保存なら転置
                if (nrow(mat) != T_expected && ncol(mat) == T_expected) {
                        mat <- t(mat)
                        message(sprintf("[BSCount] transposed to %dx%d to match T=%d.", nrow(mat), ncol(mat), T_expected))
                }
                if (nrow(mat) != T_expected) {
                        stop(sprintf("Counts matrix shape %dx%d; expected rows=%d.", nrow(mat), ncol(mat), T_expected))
                }
                storage.mode(mat) <- "double"
                mat
        }
        
        .counts_to_lambda_star <- function(counts_mat, delta_t = 1){
                if (length(delta_t) == 1) return(counts_mat / delta_t)
                sweep(counts_mat, 1, delta_t, "/")
        }
        
        # ---------- CI/帯 ----------
        .make_tx <- function(name){
                switch(name,
                       "log" = list(f = function(x) log(pmax(x, 1e-12)),   inv = function(y) exp(y)),
                       "sqrt"= list(f = function(x) sqrt(pmax(x, 0)),      inv = function(y) y^2),
                       "Anscombe" = list(f = function(x) 2*sqrt(pmax(x,0)+3/8),
                                         inv = function(y) pmax((y/2)^2 - 3/8, 0)),
                       "Fisz" = list(f = function(x) 2*sqrt(pmax(x,0)+1/4),
                                     inv = function(y) pmax((y/2)^2 - 1/4, 0)),
                       stop("Unknown transform: ", name)
                )
        }
        .ci_percentile_transform <- function(lambda_hat_star, alpha=0.05, transform="Anscombe"){
                tr <- .make_tx(transform); Y <- tr$f(lambda_hat_star)
                lo <- apply(Y, 1, quantile, probs = alpha/2,     type = 6, na.rm = TRUE)
                hi <- apply(Y, 1, quantile, probs = 1 - alpha/2, type = 6, na.rm = TRUE)
                cbind(lo = tr$inv(lo), hi = tr$inv(hi))
        }
        .ci_studentized_quick <- function(lambda_hat, lambda_hat_star, alpha=0.05){
                s_hat <- apply(lambda_hat_star, 1, sd, na.rm = TRUE); s_hat <- pmax(s_hat, 1e-12)
                Tmat  <- sweep(lambda_hat_star, 1, lambda_hat, "-")
                Tmat  <- sweep(Tmat, 1, s_hat, "/")
                qlo <- apply(Tmat, 1, quantile, probs = 1 - alpha/2, type = 6, na.rm = TRUE)
                qhi <- apply(Tmat, 1, quantile, probs = alpha/2,     type = 6, na.rm = TRUE)
                lo <- lambda_hat - qlo * s_hat; hi <- lambda_hat - qhi * s_hat
                cbind(lo = pmax(lo, 0), hi = pmax(hi, 0))
        }
        .band_supnorm <- function(lambda_hat, lambda_hat_star, alpha=0.05, studentized=FALSE){
                if (!studentized) {
                        Tb <- apply(abs(sweep(lambda_hat_star, 1, lambda_hat, "-")), 2, max, na.rm = TRUE)
                        c  <- as.numeric(quantile(Tb, probs = 1 - alpha, type = 6, na.rm = TRUE))
                        lo <- lambda_hat - c; hi <- lambda_hat + c
                } else {
                        s_hat <- apply(lambda_hat_star, 1, sd, na.rm = TRUE); s_hat <- pmax(s_hat, 1e-12)
                        Z  <- sweep(sweep(lambda_hat_star, 1, lambda_hat, "-"), 1, s_hat, "/")
                        Tb <- apply(abs(Z), 2, max, na.rm = TRUE)
                        c  <- as.numeric(quantile(Tb, probs = 1 - alpha, type = 6, na.rm = TRUE))
                        lo <- lambda_hat - c * s_hat; hi <- lambda_hat + c * s_hat
                }
                cbind(lo = pmax(lo, 0), hi = pmax(hi, 0))
        }
        
        # ---------- 出力 ----------
        .outdir_ie <- if (is_ndt) .p("IE-output", method_dir) else .p("IE-output", method_dir, variant_eff)
        
        .write_ie_each <- function(ci_mat, tag){
                dir.create(.outdir_ie, recursive = TRUE, showWarnings = FALSE)
                df <- as.data.frame(ci_mat); colnames(df)[1:2] <- c("lo","hi")
                out <- .p(.outdir_ie,
                          sprintf("DS%s-IEEachPnt-%s-%s-%s-%s-%s(J=%s).csv",
                                  DS, tag, cos_tag, kcode_eff, thrMode, thrName, J))
                write.csv(df, out, row.names = FALSE); message("[IEEachPnt] wrote: ", out)
        }
        .write_ie_band <- function(band_mat, tag){
                dir.create(.outdir_ie, recursive = TRUE, showWarnings = FALSE)
                df <- as.data.frame(band_mat); colnames(df)[1:2] <- c("lo","hi")
                out <- .p(.outdir_ie,
                          sprintf("DS%s-IEBand-%s-%s-%s-%s-%s(J=%s).csv",
                                  DS, tag, cos_tag, kcode_eff, thrMode, thrName, J))
                write.csv(df, out, row.names = FALSE); message("[IEBand] wrote: ", out)
        }
        
        # ---------- 実行 ----------
        lambda_hat <- .read_lambda_hat()
        Tlen <- length(lambda_hat)
        
        counts_mat <- .read_BSCount_strict()
        if (nrow(counts_mat) != Tlen && ncol(counts_mat) == Tlen) {
                counts_mat <- t(counts_mat)
                message(sprintf("[BSCount] transposed to %dx%d to match T=%d.", nrow(counts_mat), ncol(counts_mat), Tlen))
        }
        if (nrow(counts_mat) != Tlen) {
                stop(sprintf("Counts matrix shape %dx%d; expected rows=%d.", nrow(counts_mat), ncol(counts_mat), Tlen))
        }
        lambda_hat_star <- .counts_to_lambda_star(counts_mat, delta_t = delta_t)
        
        # 変換＋パーセンタイル
        if (length(transforms) > 0) {
                for (tr in transforms) {
                        ci <- .ci_percentile_transform(lambda_hat_star, alpha = alpha, transform = tr)
                        tag <- switch(tr,
                                      "Anscombe" = "PercAns",
                                      "log"      = "PercLog",
                                      "sqrt"     = "PercSqrt",
                                      "Fisz"     = "PercFisz",
                                      paste0("Perc", tr))
                        .write_ie_each(ci, tag)
                }
        }
        
        # 準t
        if (isTRUE(emit_tquick)) {
                ci_tq <- .ci_studentized_quick(lambda_hat, lambda_hat_star, alpha = alpha)
                .write_ie_each(ci_tq, "tQuick")
        }
        
        # 同時帯
        if (isTRUE(emit_band_sup)) {
                band <- .band_supnorm(lambda_hat, lambda_hat_star, alpha = alpha, studentized = FALSE)
                .write_ie_band(band, "BandSup")
        }
        if (isTRUE(emit_band_supt)) {
                band_t <- .band_supnorm(lambda_hat, lambda_hat_star, alpha = alpha, studentized = TRUE)
                .write_ie_band(band_t, "BandSupT")
        }
        
        invisible(list(
                method_dir = method_dir, variant = if (is_ndt) NULL else variant_eff,
                DS = DS, thrMode = thrMode, thrName = thrName, J = J,
                kcode = kcode_eff, bscount_tag = bscount_tag,
                T = Tlen, B = ncol(lambda_hat_star),
                emitted = c(if(length(transforms)) paste0("Perc:", transforms),
                            if(emit_tquick) "tQuick",
                            if(emit_band_sup) "BandSup",
                            if(emit_band_supt) "BandSupT")
        ))
}
