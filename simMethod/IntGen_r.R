Intsimulate= function(estimate,roopj,estimationMethod,thresholdRule,thresholdValue,directory_path)
{
	#シミュレーション関数の定義：生成区間 [0, t_max]
	get_nhpp_realization <- function(lambda,t_max){
		t <- 0				# シミュレーション開始時刻
		Lambda <- function(tupper){
			#.Machine$double.eps^.05=0.1649385
			integrate(f=lambda, lower=0, upper=tupper, subdivisions=200, rel.tol=.Machine$double.eps^.05)$value
		}
		Ft <- function(z) 1 - exp(-(Lambda(t+z) - Lambda(t)))
		Ft_inv <- function(u){
			a <- 0
			b <- t_max+2	# 二分探索の範囲を広げるのは，Ft の値が 1 に近づく範囲を確保するため．
			eps <- 1e-4		# 精度設定．0.0001を0とみなしている．
			while(abs(a-b)>eps){
				mid <- (a + b) / 2			# 中心点を計算
				Ft_mid <- Ft(mid)			# Ft(mid)の値を計算
				
				if (Ft_mid <= u) {
					a <- mid				# Ft(mid)がu以下なら，aをmidに更新
				} else {
					b <- mid				# Ft(mid)がuより大きいなら，bをmidに更新
				}
			}
			return(0.5*(a+b))
		}
		
		max_events <- 1000	# 最大イベント数
		X <- numeric(max_events)
		n_events <- 0		# 実際のイベント数
		
		while(t<=t_max && n_events<max_events){
			dt <- Ft_inv(runif(1))			# 次の時間間隔を計算
			if(dt > t_max){ dt <- t_max }	# 時間間隔が異常に大きい場合は t_max に修
			t <- t + dt
			n_events <- n_events + 1
			X[n_events] <- round(t,3)		# 精度設定eps<-1e-4としたため精度は小数点以下4桁までしか保証されない．
		}
		X <- X[1:n_events]	# 不要なゼロを削除．
		return(X)
	}

	#シミュレーション終了時刻の設定：推定値ベクトルの長さに等しい
	t_n <- length(estimate)

	#強度関数の補間：hat λ_i (i=1, 2, ..., t_n) ⇒ hat λ(t) (t in [0, t_n])
	timeindex <- 1:t_n
	x <- c(0,timeindex)
	y <- c(0,estimate)
	lambda <- approxfun(x, y, method="linear", f=1, rule=2)		# 線形補間(一次スプライン補間)
	#rule=1の場合，x>t_nにおけるlambda(x)はNAなので，積分エラーが起こります．
	#rule=2の場合，x>t_nにおけるlambda(x)は全てlambda(t_n)とされるので，積分エラーが起こらない．
	
	#シミュレーションの実行回数を設定
	m <- 10000
	
	#格納用変数
	arrivals_list <- list()					# 擬似時刻データt_samples × シミュレーション実行回数を格納するリスト
	counts <- matrix(0, nrow=m, ncol=t_n)	# arrivals_listを変換した擬似個数データを格納する行列(m × t_n)
	
	# シミュレーション実行：時刻データと個数データの生成
	for (i in 1:m){
		t_samples <- get_nhpp_realization(lambda, t_n)
		
		# t_samplesの中でt_nより大きい値があれば削除
		if (length(t_samples) > 0 && t_samples[length(t_samples)] > t_n){
			t_samples <- t_samples[-length(t_samples)]
		}
		
		arrivals_list[[i]] <- t_samples		# 擬似時刻データを格納
		c_samples <- hist(arrivals_list[[i]], plot=F,breaks=rep(0:t_n, by=1))$counts
		counts[i,] <- c_samples				# 擬似個数データを格納
	}
	
	#リストarrivals_listを行列に変換
	max_length <- max(sapply(arrivals_list, length))			# arrivals_listの要素の最大長を取得
	arrivals <- matrix(NA, nrow=m, ncol=max_length)				# 最大長に基づいた行列の初期化
	for (i in 1:m) {
		arrivals[i, 1:length(arrivals_list[[i]])] <- arrivals_list[[i]]
	}
	
	#擬似時刻データ出力
	write.table(arrivals, file = "arrivals.dat", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
	#擬似個数データ出力
	write.table(counts, file = "counts.dat", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}