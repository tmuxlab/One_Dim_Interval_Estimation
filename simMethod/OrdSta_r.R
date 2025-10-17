Ordsimulate= function(estimate,roopj,estimationMethod,thresholdRule,thresholdValue,directory_path)
{
  #シミュレーション終了時刻の設定：推定値ベクトルの長さに等しい
  t_n <- length(estimate)
  t_max <- t_n
  
  #強度関数の補間：hat λ_i (i=1, 2, ..., t_n) ⇒ hat λ(t) (t in [0, t_n])
  timeindex <- 1:t_n
  x <- c(0,timeindex)
  y <- c(0,estimate)
  
  ###ここ，元々rule=1でやってましたが，rule=2に統一した方がいいですか． ###
  
  lambda <- approxfun(x, y, method="linear", f=1, rule=2)#線形補間(一次スプライン補間)
  #rule=1の場合，x>t_nにおけるlambda(x)はNAなので，積分エラーが起こります．
  #rule=2の場合，x>t_nにおけるlambda(x)は全てlambda(t_n)とされるので，積分エラーが起こらない．
  Lambda <- function(tupper) { #.Machine$double.eps^.05=0.1649385
          integrate(f=lambda, lower=0, upper=tupper, subdivisions=200, rel.tol=.Machine$double.eps^.05)$value
  }
  total_mass <- Lambda(t_n)
  # 安全策：質量が非正ならゼロ到着（必要に応じて分岐を追加）
  if (!is.finite(total_mass) || total_mass <= 0) {
    m <- 10000
    counts <- matrix(0, nrow = m, ncol = t_n)
    arrivals <- matrix(NA, nrow = m, ncol = 1)
    write.table(arrivals, paste0(dirname(rstudioapi::getSourceEditorContext()$path), directory_path,
                                 "BS-OrdSta-", estimationMethod, "-", thresholdRule, "-", thresholdValue, "(J=", roopj, ").dat"),
                append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE)
    write.table(counts,   paste0(dirname(rstudioapi::getSourceEditorContext()$path), directory_path,
                                 "BSCount-OrdSta-", estimationMethod, "-", thresholdRule, "-", thresholdValue, "(J=", roopj, ").dat"),
                append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE)
    return(invisible(NULL))
  }
  
  Ft <- function(z) Lambda(z)/total_mass #DS2での実行不可問題対処考察中
  #Ft <- function(z) sapply(z,Lambda(z))/sapply(t_max,Lambda)
  
  Ft_inv <- function(u){
    a <- 0
    b <- t_max+2              # 二分探索の範囲を広げるのは，Ft の値が 1 に近づく範囲を確保するため．
    eps <- 1e-4               # 精度設定．0.0001を0とみなしている．
    while(abs(a-b)>eps){
      mid <- (a + b)/2  # 中心点を計算
      Ft_mid <- Ft(mid) # Ft(mid)の値を計算
      if(Ft_mid<=u){
        a <- mid  # Ft(mid)がu以下なら，aをmidに更新
      }else{
        b <- mid  # Ft(mid)がuより大きいなら，bをmidに更新
      }
    }
    return(0.5*(a+b))
  }
  
  #シミュレーション関数の定義：生成区間 [0, t_max]
  get_nhpp_realization <- function(Ft_inv,total_mass){
    n <- rpois(1, total_mass)
    if (n <= 0) return(numeric(0))       # ← 安全化
    X <- sapply(1:n, function(v) Ft_inv(runif(1)))
    return(sort(X))
  }
  
  #シミュレーションの実行回数を設定
  m <- 10000
  # 格納用変数
  arrivals_list <- list()                # 擬似時刻データt_samples × シミュレーション実行回数を格納するリスト
  counts <- matrix(0, nrow=m, ncol=t_n)  # arrivals_listを変換した擬似個数データを格納する行列(m × t_n)
  
  # シミュレーション実行：時刻データと個数データの生成
  for (i in 1:m){
          t_samples <- get_nhpp_realization(Ft_inv,total_mass)
          # t_samplesの中でt_nより大きい値を削除する
          if (length(t_samples) > 0 && t_samples[length(t_samples)] > t_n) {
                  t_samples <- t_samples[-length(t_samples)]
          }
          arrivals_list[[i]]<- t_samples       # 擬似時刻データを格納
          c_samples <- hist(arrivals_list[[i]], plot=F,breaks = 0:t_n)$counts
          counts[i,] <- c_samples              # 擬似個数データを格納
  }
  
  #リストarrivals_listを行列に変換
  max_length <- max(sapply(arrivals_list, length))        # arrivals_listの要素の最大長を取得
  arrivals <- matrix(NA, nrow=m, ncol=max_length)         # 最大長に基づいた行列の初期化
  for (i in 1:m) {
          arrivals[i, 1:length(arrivals_list[[i]])] <- arrivals_list[[i]]
  }
  
  #時間データ出力
  write.table(arrivals,paste0(dirname(rstudioapi::getSourceEditorContext()$path),directory_path,"BS-OrdSta-",estimationMethod,"-",thresholdRule,"-",thresholdValue,"(J=",roopj,").dat"), append=F, quote=F, col.names=F, row.names=F)#出力ファイル名
  #個数データ出力
  write.table(counts,paste0(dirname(rstudioapi::getSourceEditorContext()$path),directory_path,"BSCount-OrdSta-",estimationMethod,"-",thresholdRule,"-",thresholdValue,"(J=",roopj,").dat"), append=F, quote=F, col.names=F, row.names=F)#出力ファイル名
}
