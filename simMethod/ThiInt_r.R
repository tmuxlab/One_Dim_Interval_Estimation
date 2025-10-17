ThiIntsimulate= function(estimate,roopj,estimationMethod,thresholdRule,thresholdValue,directory_path)
{
  #シミュレーション関数の定義：生成区間 [0, t_max]
  get_nhpp_realization <- function(t_max,lambda_star_value){
    t <- 0                     # シミュレーション開始時刻
    X <- numeric(0)            # 到着時刻を格納するベクトル
    while(t <= t_max){
        # 積分法（inverse transform method）によって次の到着間隔を生成
        # Ft_inv は累積強度関数 Λ(t) の逆関数（またはその近似）であり、
        # 一様乱数から非定常ポアソン過程の inter-arrival time をサンプリングする
        t <- t + Ft_inv(runif(1))
        
        # 到着時刻が t_max を超えたらループ終了（lambda(t) の定義外を防ぐ）
        if(t > t_max) break
        
        # 棄却法（thinning method）の要素
        # Ft_inv が近似である可能性があるため、補正として到着時刻 t を
        # 確率 lambda(t)/lambda_star() で受け入れる
        # これにより、最終的に λ(t) に従った非定常ポアソン過程が得られる
        if(runif(1) < lambda(t) / lambda_star_value) {
              X <- c(X, t)           # 受け入れた到着時刻を記録
        }
    }
    return(X)
  }
        
  #シミュレーション終了時刻の設定：推定値ベクトルの長さに等しい
  t_n <- length(estimate)
  #強度関数の補間：hat λ_i (i=1, 2, ..., t_n) ⇒ hat λ(t) (t in [0, t_n])
  timeindex <- 1:t_n
  x <- c(0,timeindex)
  y <- c(0,estimate)
  
  lambda <- approxfun(x, y, method="linear", f=1, rule=2)#線形補間(一次スプライン補間)
  #rule=1の場合，x>t_nにおけるlambda(x)はNAなので，積分エラーが起こります．
  #rule=2の場合，x>t_nにおけるlambda(x)は全てlambda(t_n)とされるので，積分エラーが起こらない．
  lambda_star_value <- max(sapply(seq(0, t_n, length.out = 1000), lambda))
  Ft_inv <- function(u) {-log(1 - u) / lambda_star_value}

  #シミュレーションの実行回数を設定
  m <- 10000
  # 格納用変数
  arrivals_list <- list()                # 擬似時刻データt_samples × シミュレーション実行回数を格納するリスト
  counts <- matrix(0, nrow=m, ncol=t_n)  # arrivals_listを変換した擬似個数データを格納する行列(m × t_n)
  
  # シミュレーション実行：時刻データと個数データの生成
  for (i in 1:m){
          t_samples <- get_nhpp_realization(t_n,lambda_star_value)
          # t_samplesの中でt_nより大きい値を削除する
          if (length(t_samples) > 0 && t_samples[length(t_samples)] > t_n) {
                  t_samples <- t_samples[-length(t_samples)]
          }
          arrivals_list[[i]]<- t_samples       # 擬似時刻データを格納
          c_samples <- hist(arrivals_list[[i]], plot=F,breaks=0:t_n)$counts
          counts[i,] <- c_samples              # 擬似個数データを格納
  }
  
  #リストarrivals_listを行列に変換
  max_length <- max(sapply(arrivals_list, length))        # arrivals_listの要素の最大長を取得
  arrivals <- matrix(NA, nrow=m, ncol=max_length)         # 最大長に基づいた行列の初期化
  for (i in 1:m) {
          arrivals[i, 1:length(arrivals_list[[i]])] <- arrivals_list[[i]]
  }
  
  #時間データ出力
  write.table(arrivals,paste0(dirname(rstudioapi::getSourceEditorContext()$path),directory_path,"BS-ThiInt-",estimationMethod,"-",thresholdRule,"-",thresholdValue,"(J=",roopj,").dat"), append=F, quote=F, col.names=F, row.names=F)#出力ファイル名
  #個数データ出力
  write.table(counts,paste0(dirname(rstudioapi::getSourceEditorContext()$path),directory_path,"BSCount-ThiInt-",estimationMethod,"-",thresholdRule,"-",thresholdValue,"(J=",roopj,").dat"), append=F, quote=F, col.names=F, row.names=F)#出力ファイル名
}
