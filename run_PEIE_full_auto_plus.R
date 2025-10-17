# ============================================================
# フルオート拡張：生成 → 可視化(帯2種) → 検定(PIT/カバレッジ) → PITヒスト&QQ
# 依存：run_PEIE_addons_generic() が読み込み済みであること
# ============================================================

# ---- 単条件（DS × method_dir/variant × COS1件 × J1件）を処理する本体 ----
run_PEIE_full_auto_ex <- function(
                DS, thrMode, thrName, J,
                method_dir, variant = NULL,             # NDT_WSE は variant=NULL、DT系は A1/A2/A3 or B1/B2
                cos_tag = "COS",
                transforms = c("Anscombe","log","sqrt"),
                alpha = 0.05, delta_t = 1,
                emit_tquick = TRUE, emit_band_sup = TRUE, emit_band_supt = TRUE
){
        stopifnot(exists("run_PEIE_addons_generic"))
        # 1) 追加CI/帯を生成（CSV出力）
        run_PEIE_addons_generic(
                DS=DS, thrMode=thrMode, thrName=thrName, J=J,
                method_dir=method_dir, variant=variant, cos_tag=cos_tag,
                transforms=transforms, alpha=alpha, delta_t=delta_t,
                emit_tquick=emit_tquick, emit_band_sup=emit_band_sup, emit_band_supt=emit_band_supt
        )
        
        # 2) 可視化＆検定
        is_ndt <- identical(method_dir, "NDT_WSE")
        kcode_eff <- if (is_ndt) "kH" else paste0("k", variant)
        bscount_tag <- if (is_ndt) "H" else variant
        
        .p <- function(...) { parts <- list(...); parts <- parts[!vapply(parts, function(x) is.null(x) || (is.character(x)&&x==""), logical(1))]; do.call(file.path, parts) }
        stage_dir <- function(stage){ if (is_ndt) .p(stage, method_dir) else .p(stage, method_dir, variant) }
        
        # --- 読み込み（λ̂, CI, 帯, BSCount→λ*） ---
        read_lambda_hat <- function(){
                f <- .p(stage_dir("output"), sprintf("DS%s-%s-%s-%s(J=%s).txt", DS, kcode_eff, thrMode, thrName, J))
                df <- try(read.csv(f, header=TRUE, check.names=FALSE, strip.white=TRUE, stringsAsFactors=FALSE), silent=TRUE)
                if (inherits(df,"try-error")) df <- try(read.table(f, header=TRUE, check.names=FALSE, strip.white=TRUE, stringsAsFactors=FALSE), silent=TRUE)
                if (inherits(df,"try-error")) df <- read.table(f, header=FALSE, strip.white=TRUE)
                as.numeric(df[[1]])
        }
        read_ci <- function(tag){
                f <- .p(stage_dir("IE-output"), sprintf("DS%s-IEEachPnt-%s-%s-%s-%s-%s(J=%s).csv", DS, tag, cos_tag, kcode_eff, thrMode, thrName, J))
                if (!file.exists(f)) return(NULL)
                df <- read.csv(f, header=TRUE, check.names=FALSE); if (!all(c("lo","hi") %in% names(df))) return(NULL)
                df
        }
        read_band <- function(tag){
                f <- .p(stage_dir("IE-output"), sprintf("DS%s-IEBand-%s-%s-%s-%s-%s(J=%s).csv", DS, tag, cos_tag, kcode_eff, thrMode, thrName, J))
                if (!file.exists(f)) return(NULL)
                df <- read.csv(f, header=TRUE, check.names=FALSE); if (!all(c("lo","hi") %in% names(df))) return(NULL)
                df
        }
        read_BSCount <- function(T_expected){
                f <- .p(stage_dir("BS-output"), sprintf("DS%s-BSCount-%s-%s-%s-%s(J=%s).dat", DS, cos_tag, bscount_tag, thrMode, thrName, J))
                df <- try(read.table(f, header=TRUE, check.names=FALSE, strip.white=TRUE, stringsAsFactors=FALSE), silent=TRUE)
                if (inherits(df,"try-error")) df <- read.table(f, header=FALSE, check.names=FALSE, strip.white=TRUE)
                df[] <- lapply(df, function(col) suppressWarnings(as.numeric(col)))
                mat <- as.matrix(df); storage.mode(mat) <- "double"
                if (nrow(mat)!=T_expected && ncol(mat)==T_expected) { mat <- t(mat); message(sprintf("[BSCount %s] transposed to %dx%d (T=%d).", cos_tag, nrow(mat), ncol(mat), T_expected)) }
                if (nrow(mat)!=T_expected) stop(sprintf("BSCount(%s) shape %dx%d; expected rows=%d.", cos_tag, nrow(mat), ncol(mat), T_expected))
                mat
        }
        counts_to_lambda_star <- function(counts_mat, delta_t=1){
                if (length(delta_t)==1) counts_mat/delta_t else sweep(counts_mat, 1, delta_t, "/")
        }
        
        lambda_hat <- read_lambda_hat(); Tlen <- length(lambda_hat)
        ciA <- read_ci("PercAns"); ciL <- read_ci("PercLog"); ciS <- read_ci("PercSqrt"); ciT <- if (emit_tquick) read_ci("tQuick") else NULL
        band  <- if (emit_band_sup)  read_band("BandSup")  else NULL
        bandT <- if (emit_band_supt) read_band("BandSupT") else NULL
        
        counts <- read_BSCount(T_expected=Tlen)
        lambda_star <- counts_to_lambda_star(counts, delta_t = delta_t)
        
        # --- 可視化（BandSup & BandSupT を帯で重ね、CI線とλ̂を表示） ---
        viz_dir <- stage_dir("IE-output")
        pngfile <- .p(viz_dir, sprintf("DS%s-viz-%s-%s-%s-%s(J=%s).png", DS, cos_tag, kcode_eff, thrMode, thrName, J))
        png(pngfile, width=1400, height=900, res=150)
        par(mar=c(4,5,2,1))
        plot(lambda_hat, type="l", lwd=2, xlab="t (index/day)", ylab="lambda", main=sprintf("%s: CIs + Bands", cos_tag))
        if (!is.null(band)) {
                polygon(c(seq_len(Tlen), rev(seq_len(Tlen))), c(band$lo, rev(band$hi)), border=NA, col=rgb(0.6,0.6,0.6,0.35))
        }
        if (!is.null(bandT)) {
                polygon(c(seq_len(Tlen), rev(seq_len(Tlen))), c(bandT$lo, rev(bandT$hi)), border="#555555", col=rgb(0.3,0.3,0.3,0.18))
        }
        if (!is.null(ciA)) { lines(ciA$lo, lty=2); lines(ciA$hi, lty=2) }
        if (!is.null(ciL)) { lines(ciL$lo, lty=4); lines(ciL$hi, lty=4) }
        if (!is.null(ciS)) { lines(ciS$lo, lty=5); lines(ciS$hi, lty=5) }
        if (!is.null(ciT)) { lines(ciT$lo, col=4, lty=3); lines(ciT$hi, col=4, lty=3) }
        lines(lambda_hat, col=2, lwd=2)
        legend("topright",
               legend=c("lambda_hat","PercAns","PercLog","PercSqrt","tQuick","BandSup","BandSupT"),
               lty=c(1,2,4,5,3,NA,NA), lwd=c(2,1,1,1,1,NA,NA),
               pch=c(NA,NA,NA,NA,NA,15,15), col=c(2,1,1,1,4,gray(0.6),gray(0.3)), bty="n", pt.cex=2)
        dev.off(); message("Saved plot: ", pngfile)
        
        # --- 検定：PIT(KS) と Coverage（二項） ---
        make_tx <- function(name){
                switch(name,
                       "log"      = list(f=function(x) log(pmax(x,1e-12))),
                       "sqrt"     = list(f=function(x) sqrt(pmax(x,0))),
                       "Anscombe" = list(f=function(x) 2*sqrt(pmax(x,0)+3/8))
                )
        }
        pit_from_boot <- function(theta_hat, theta_star, transform=c("none","Anscombe","log","sqrt")){
                transform <- match.arg(transform); th <- theta_hat; ts <- theta_star
                if (transform!="none"){ tr <- make_tx(transform); th <- tr$f(th); ts <- tr$f(ts) }
                pits <- numeric(length(th)); for (t in seq_along(th)){ x <- ts[t,]; r <- sum(x <= th[t], na.rm=TRUE); pits[t] <- r/(sum(!is.na(x))+1) }; pits
        }
        coverage_test <- function(lambda_hat, ci_df, alpha=0.05){
                inside <- (lambda_hat >= ci_df$lo) & (lambda_hat <= ci_df$hi)
                S <- sum(inside, na.rm=TRUE); Tlen <- length(lambda_hat); p <- binom.test(S, Tlen, p=1-alpha)$p.value
                list(S=S, T=Tlen, prop=S/Tlen, p=p)
        }
        
        # PIT: PercAns / tQuick
        pit_ans <- if (!is.null(ciA)) pit_from_boot(lambda_hat, lambda_star, transform="Anscombe") else NULL
        pit_tq  <- if (!is.null(ciT)) {
                s_hat <- apply(lambda_star, 1, sd, na.rm=TRUE); s_hat <- pmax(s_hat, 1e-12)
                Z <- sweep(sweep(lambda_star, 1, lambda_hat, "-"), 1, s_hat, "/")
                pits <- numeric(Tlen); for (t in 1:Tlen){ z <- Z[t,]; r <- sum(z <= 0, na.rm=TRUE); pits[t] <- r/(sum(!is.na(z))+1) }; pits
        } else NULL
        
        ks_ans <- if (!is.null(pit_ans)) suppressWarnings(ks.test(pit_ans, "punif", 0, 1)) else NULL
        ks_tq  <- if (!is.null(pit_tq))  suppressWarnings(ks.test(pit_tq,  "punif", 0, 1)) else NULL
        
        cov_ans <- if (!is.null(ciA)) coverage_test(lambda_hat, ciA, alpha) else NULL
        cov_tq  <- if (!is.null(ciT)) coverage_test(lambda_hat, ciT, alpha) else NULL
        
        # PITヒスト＆QQ（PercAns / tQuick）
        pit_plot <- function(pits, label){
                if (is.null(pits)) return(NULL)
                png(.p(viz_dir, sprintf("DS%s-viz-PIT-%s-%s-%s-%s-%s(J=%s).png", DS, label, cos_tag, kcode_eff, thrMode, thrName, J)),
                    width=1400, height=600, res=150)
                par(mfrow=c(1,2), mar=c(4,4,2,1))
                hist(pits, breaks=10, main=sprintf("PIT Hist (%s, %s)", label, cos_tag), xlab="PIT", col="grey80", border="white")
                abline(h=length(pits)/10, lty=2, col="grey40")
                u <- (1:length(pits)-0.5)/length(pits)
                plot(sort(u), sort(pits), xlab="Theoretical U(0,1)", ylab="Empirical PIT", main=sprintf("PIT QQ (%s, %s)", label, cos_tag), pch=19, cex=.7)
                abline(0,1,col="red")
                dev.off()
        }
        pit_plot(pit_ans, "PercAns"); pit_plot(pit_tq, "tQuick")
        
        # サマリー1件を返す
        data.frame(
                DS=DS, Method=method_dir, Variant=if (is_ndt) "H" else variant,
                COS=cos_tag, J=J, thrMode=thrMode, thrName=thrName,
                KS_PercAns_D = if (!is.null(ks_ans)) unname(ks_ans$statistic) else NA_real_,
                KS_PercAns_p = if (!is.null(ks_ans)) ks_ans$p.value else NA_real_,
                KS_tQuick_D  = if (!is.null(ks_tq))  unname(ks_tq$statistic) else NA_real_,
                KS_tQuick_p  = if (!is.null(ks_tq))  ks_tq$p.value else NA_real_,
                Cov_PercAns_inside = if (!is.null(cov_ans)) cov_ans$S else NA_integer_,
                Cov_PercAns_prop   = if (!is.null(cov_ans)) cov_ans$prop else NA_real_,
                Cov_PercAns_p      = if (!is.null(cov_ans)) cov_ans$p else NA_real_,
                Cov_tQuick_inside  = if (!is.null(cov_tq))  cov_tq$S else NA_integer_,
                Cov_tQuick_prop    = if (!is.null(cov_tq))  cov_tq$prop else NA_real_,
                Cov_tQuick_p       = if (!is.null(cov_tq))  cov_tq$p else NA_real_,
                stringsAsFactors = FALSE
        )
}

# --- 追加: 失敗時にNAで埋めた1行を返すテンプレ ---
.make_na_row <- function(DS, thrMode, tn, J, md, vr, ct){
        data.frame(
                DS = DS, Method = md, Variant = if (identical(md,"NDT_WSE")) "H" else vr,
                COS = ct, J = J, thrMode = thrMode, thrName = tn,
                KS_PercAns_D = NA_real_, KS_PercAns_p = NA_real_,
                KS_tQuick_D  = NA_real_, KS_tQuick_p  = NA_real_,
                Cov_PercAns_inside = NA_integer_, Cov_PercAns_prop = NA_real_, Cov_PercAns_p = NA_real_,
                Cov_tQuick_inside  = NA_integer_, Cov_tQuick_prop  = NA_real_, Cov_tQuick_p  = NA_real_,
                stringsAsFactors = FALSE
        )
}

# ---- 総当たりランナー（COS × J × 体制） ----
run_PEIE_all <- function(
                DS, thrMode,
                J_vec,                                          # 例: 2:5
                tasks = list(                                   # 体制を定義
                        list(method_dir="NDT_WSE",    variant=NULL, thrName="ldt"), # NDT_H（ldt専用）
                        list(method_dir="DT_Ans_WSE", variant="A3", thrName="ut"),  # 例: A3
                        list(method_dir="DT_Bar_WSE", variant="B2", thrName="ut")   # 例: B2
                ),
                cos_tags = c("COS","GS-","GS+","OrdSta","ThiArr","ThiInt"),
                transforms = c("Anscombe","log","sqrt"),
                alpha = 0.05, delta_t = 1
){
        out_all <- list()
        for (tsk in tasks){
                md <- tsk$method_dir; vr <- tsk$variant; tn <- tsk$thrName
                for (J in J_vec){
                        for (ct in cos_tags){
                                # --- ここを tryCatch でガード ---
                                res <- tryCatch(
                                        run_PEIE_full_auto_ex(
                                                DS=DS, thrMode=thrMode, thrName=tn, J=J,
                                                method_dir=md, variant=vr, cos_tag=ct,
                                                transforms=transforms, alpha=alpha, delta_t=delta_t,
                                                emit_tquick=TRUE, emit_band_sup=TRUE, emit_band_supt=TRUE
                                        ),
                                        error = function(e){
                                                # 見つからない等の理由で失敗 -> コメント出してNA行を返す
                                                message(sprintf(
                                                        "[SKIP] DS=%s J=%s %s/%s COS=%s : %s",
                                                        DS, J, md, if (is.null(vr)) "H" else vr, ct, conditionMessage(e)
                                                ))
                                                # ヒートマップ等の形を保つため、NAを入れた行を作る
                                                .make_na_row(DS, thrMode, tn, J, md, vr, ct)
                                        }
                                )
                                # run_PEIE_full_auto_ex は成功時 data.frame を返すのでそのまま積む
                                out_all[[length(out_all)+1]] <- res
                        }
                }
        }
        
        # 失敗分を含めて結合（NAはそのまま）
        summary_df <- do.call(rbind, out_all)
        
        # 以降は元のまま（体制ごとのCSV出力）
        write_summary <- function(md, vr, df){
                base_dir <- if (identical(md,"NDT_WSE")) file.path("IE-output", md) else file.path("IE-output", md, vr)
                dir.create(base_dir, recursive=TRUE, showWarnings = FALSE)
                fn <- file.path(base_dir, sprintf("DS%s-GRID-SUMMARY-%s(thrMode=%s).csv", unique(df$DS), md, unique(df$thrMode)))
                write.csv(df[df$Method==md & (if (md=="NDT_WSE") TRUE else df$Variant==vr), ], fn, row.names = FALSE)
                message("Saved grid summary: ", fn)
        }
        for (tsk in tasks){
                write_summary(tsk$method_dir, if (is.null(tsk$variant)) "H" else tsk$variant, summary_df)
        }
        invisible(summary_df)
}
