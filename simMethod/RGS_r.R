RGSsimulate= function(estimate,roopj,estimationMethod,thresholdRule,thresholdValue,directory_path)
{
  tn <- length(estimate)
  #入力区間 [1,tn]
  timeindex <- 1:tn
  #補間 ^λ(t)
  #補間区間 [0,tn]
  x <- c(0,timeindex)
  y <- c(0,estimate)
  #格納行列
  row <- 10000
  col <- 500
  # arrivals 格納行列 : row * col
  arrivals <- matrix(0, nrow=row, ncol=col)
  # counts 格納行列 ： row * tn
  counts <- matrix(0, nrow=row, ncol=tn)
  
  # y > 0 のデータに限定
  valid_idx <- which(y > 0)
  x_valid <- x[valid_idx]
  y_valid <- y[valid_idx]
  
  # ログ変換
  logy <- log(y_valid)
  
  # 線形回帰で係数を推定
  lm_fit <- lm(logy ~ x_valid + I(x_valid^2))
  
  # 初期値を取得
  start_vals <- as.list(coef(lm_fit))
  names(start_vals) <- c("a", "b", "c")
  
  # NLS に渡す
  RGS <- nls(y ~ exp(a + b * x + c * x^2),
             start = start_vals,
             algorithm = "port",
             lower = list(a = -100, b = -100, c = -100),
             upper = list(a = 100, b = 100, c = 100))
  #確認用コード
  #plot(x, y, main = "Estimate and Fitted λ(t)", col = "blue")
  #lines(x, predict(RGS), col = "red")
  #legend("topleft", legend = c("Estimate", "Fitted"), col = c("blue", "red"), lty = 1)
  
  alpha0 <- coefficients(RGS)[[1]]
  alpha1 <- coefficients(RGS)[[2]]
  alpha2 <- coefficients(RGS)[[3]]
  
  #simulation algorithm
  #生成区間 [0,tn]
  get_nhpp_realization <- function(alpha0,alpha1,alpha2){
    #step 1
    params <- step1_parameter_setting(alpha0, alpha1, alpha2, tn)
    gamma0 <- params$gamma0
    gamma1 <- params$gamma1
    a <- params$a
    b <- params$b
    #アルゴリズム１の実行
    if(gamma1<0){
      T <- algorithm1(gamma0,gamma1)
      m1 <- length(T)
      result1 <- T
    }
    if(gamma1>0){
      T <- algorithm1dash(gamma0,gamma1)
      m1 <- length(T)
      result1 <- T
    }
    #アルゴリズム２の実行
    c <- csetting(alpha0,alpha1,alpha2,b,tn)
    T <- algorithm2(alpha0,alpha1,alpha2,gamma0,gamma1,a,b,c)
    n1 <- length(T)
    result2 <- T
    #bの値によって終了判定
    if(b==tn){
      if(m1==0 && n1==0){
        return(0)
      }
      if(n1==0){
        return(result1)
      }
      if(m1==0){
        return(result2)
      }
      if(m1!=0 && n1!=0){
        #len1 <- length(result1)
        #len2 <- length(result2)
        #if(len1<len2){
        #        result1 <- c(result1,rep(0,len2-len1))
        #}
        #if(len1>len2){
        #        result2 <- c(result2,rep(0,len1-len2))
        #}
        result <- c(result1,result2)
        result <- sort(result)
        return(result)
      }
    }
    #(a,b)=(-a1/2a2,tn)として計算続行
    a <- b
    b <- tn
    #step 6
    params2 <- step6_parameter_setting(alpha0, alpha1, alpha2, tn)
    gamma0 <- params2$gamma0
    gamma1 <- params2$gamma1
    #アルゴリズム１の実行
    if(gamma1<0){
      T <- algorithm1(gamma0,gamma1)
      m2 <- length(T)
      result1dash <- T
    }
    if(gamma1>0){
      T <- algorithm1dash(gamma0,gamma1)
      m2 <- length(T)
      result1dash <- T
    }
    #アルゴリズム２の実行
    cdash <- csettingdash(alpha0,alpha1,alpha2,b,tn)
    T <- algorithm2(alpha0,alpha1,alpha2,gamma0,gamma1,a,b,cdash)
    n2 <- length(T)
    result2dash <- T
    m <- m1 + m2
    n <- n1 + n2
    result1 <- c(result1,result1dash)
    result2 <- c(result2,result2dash)
    if(m ==0 && n==0){
      return(0)
    }
    if(n==0){
      return(result1)
    }
    if(m==0){
      result2 <- sort(result2)
      return(result2)
    }
    if(m!=0 && n!=0){
      result <- c(result1,result2)
      result <- sort(result)
      return(result)
      #連結した結果(上),mergeが＋の意味(下)
      #len1 <- length(result1)
      #len2 <- length(result2)
      #if(len1<len2){
      #        result1 <- c(result1,rep(0,len2-len1))
      #}
      #if(len1>len2){
      #        result2 <- c(result2,rep(0,len1-len2))
      #}
      #return(result1+result2)
    }
    
  }
  
  step1_parameter_setting <- function(alpha0, alpha1, alpha2, t0) {
          # 初期区間 (a,b)
          a <- 0
          b <- t0
          
          # b の条件：α1·α2 < 0 かつ -α1/(2α2) < t0 のとき、b = -α1/(2α2)
          if (alpha1 * alpha2 < 0 && -alpha1 / (2 * alpha2) < t0) {
                  b <- -alpha1 / (2 * alpha2)
          }
          
          # γ0 と γ1 の設定
          if (alpha2 > 0) {
                  if (alpha1 >= 0) {
                          gamma0 <- alpha0 - alpha2 * t0^2
                          gamma1 <- alpha1 + 2 * alpha2 * t0
                  } else {
                          gamma0 <- alpha0
                          gamma1 <- alpha1
                  }
          } else if (alpha2 < 0) {
                  gamma0 <- alpha0
                  if (-alpha1 / (2 * alpha2) < t0 && alpha1 > 0) {
                          gamma1 <- alpha1 / 2
                  } else {
                          gamma1 <- alpha1 + alpha2 * t0
                  }
          } else {
                  # α2 = 0 の場合のフォールバック（通常は想定しないが安全のため）
                  gamma0 <- alpha0
                  gamma1 <- alpha1
          }
          
          # 結果をリストで返す
          return(list(gamma0 = gamma0, gamma1 = gamma1, a = a, b = b))
  }
  
  step6_parameter_setting <- function(alpha0, alpha1, alpha2, t0) {
          if (alpha1 > 0 && alpha2 < 0) {
                  gamma0 <- alpha0 + (alpha1 / 2) * t0
                  gamma1 <- (alpha1 / 2) + alpha2 * t0
          } else if (alpha1 < 0 && alpha2 > 0) {
                  gamma0 <- alpha0 - alpha2 * t0^2
                  gamma1 <- alpha1 + 2 * alpha2 * t0
          } else {
                  # それ以外の組み合わせがある場合は元の値をそのまま使う（または明示的にエラーを出してもよい）
                  gamma0 <- alpha0
                  gamma1 <- alpha1
          }
          
          return(list(gamma0 = gamma0, gamma1 = gamma1))
  }
  
  algorithm1 <- function(gamma0,gamma1){
    lambda <- exp(gamma0)
    sita <- -lambda/gamma1
    beta <- -gamma1
    i <- 1                      #ループ変数，到着時刻Tiのi.=1,2,3...m
    t <- 0                      #求めている到着時刻Ti
    m <- rpois(1,sita)          #平均Θ=-λ/α1を持つポアソン変量m
    A <- numeric(0)             #出力したいイベント時刻Tiの列
    while(i <= m){
      t <- rexp(1,rate=1)/beta/(m-i+1) + t
      if(t > tn) break          # t > tn の時，生成する時間列が区間を超えるため break
      A <- c(A,t)               #イベント時刻Tiを追加．
      i <- i + 1
    }
    return(A)
  }
  
  algorithm1dash <- function(gamma0,gamma1){
    lambda <- exp(gamma0)
    alpha <- -gamma1
    sita <- -lambda/alpha
    beta <- gamma1
    i <- 1                      #ループ変数，到着時刻Tiのi.=1,2,3...m
    t <- 0                      #求めている到着時刻Ti
    m <- rpois(1,sita)          #平均Θ=-λ/α1を持つポアソン変量m
    A <- numeric(0)             #出力したいイベント時刻Tiの列
    while(i <= m){
      t <- rexp(1,rate=1)/beta/(m-i+1) + t
      if(t > tn) break          # t > tn の時，生成する時間列が区間を超えるため break
      A <- c(A,t)               #イベント時刻Tiを追加．
      i <- i + 1
    }
    B <- A
    i <- 1
    LL <- length(A)
    while(i<=LL){
      A[i] <- tn - B[LL+1-i]
      i <- i + 1
    }
    return(A)
  }
  
  algorithm2 <- function(alpha0,alpha1,alpha2,gamma0,gamma1,a,b,c){
    k <- 0
    t <- 0
    A <- numeric(0)
    mean <- meangetting(alpha0,alpha1,alpha2,gamma0,gamma1,a,b)
    n <- rpois(1,mean)
    if(n==0){return(numeric(0))}
    while(k<n){
      V <- runif(1)
      T <- a + (b-a)*V
      U <- runif(1)
      if(U<=lambda_star(alpha0,alpha1,alpha2,gamma0,gamma1,T)/c){
        A <- c(A,T)
      }
      k <- k + 1
    }
    return(A)
  }
  
  csetting <- function(alpha0, alpha1, alpha2, b, tn){
          if(alpha1 >= 0 && alpha2 > 0){
                  c <- exp(alpha0) - exp(alpha0 - alpha2 * b^2)
          } else if(alpha1 < 0 && alpha2 > 0 && -alpha1 / (2 * alpha2) >= tn){
                  c <- exp(alpha0 + alpha1 * b + alpha2 * b^2) - exp(alpha0 + alpha1 * b)
          } else if(alpha1 < 0 && alpha2 > 0 && -alpha1 / (2 * alpha2) < tn){
                  c <- exp(alpha0 - (alpha1^2) / (4 * alpha2)) - exp(alpha0 - (alpha1^2) / (2 * alpha2))
          } else if(alpha1 <= 0 && alpha2 < 0){
                  c <- exp(alpha0) - exp(alpha0 + alpha1 * b + alpha2 * b^2)
          } else if(alpha1 > 0 && alpha2 < 0 && -alpha1 / (2 * alpha2) > tn){
                  c <- exp(alpha0 + alpha1 * b + alpha2 * b^2) - exp(alpha0)
          } else if(alpha1 > 0 && alpha2 < 0 && -alpha1 / (2 * alpha2) < tn){
                  c <- exp(alpha0 - (alpha1^2) / (4 * alpha2)) - exp(alpha0)
          } else {
                  stop("No matching case for csetting()")
          }
          return(c)
  }
  
  csettingdash <- function(alpha0, alpha1, alpha2, b){
          if(alpha1 > 0 && alpha2 < 0){
                  c <- exp(alpha0 - (alpha1^2) / (4 * alpha2)) - exp(alpha0 + alpha1 * b + alpha2 * b^2)
          } else if(alpha1 < 0 && alpha2 > 0){
                  c <- exp(alpha0 - (alpha1^2) / (4 * alpha2)) - exp(alpha0 - (alpha1^2) / (2 * alpha2) - alpha1 * b - alpha2 * b^2)
          } else {
                  stop("No matching case for csettingdash()")
          }
          return(c)
  }
  
  meangetting <- function(alpha0,alpha1,alpha2,gamma0,gamma1,a,b){
    lambda <- function(t){lambda_star(alpha0,alpha1,alpha2,gamma0,gamma1,t)}
    mean <- integrate(f=lambda,lower=a,upper=b,rel.tol=.Machine$double.eps^.05)$value
    return(mean)
  }
  
  lambda_star <- function(alpha0,alpha1,alpha2,gamma0,gamma1,T){
    result <- exp(alpha0+alpha1*T+alpha2*T**2)-exp(gamma0+gamma1*T)
    return(max(result,0))       #負の値を０に補正
  }
  
  #時間,個数データ生成
  for (i in 1:row){
    t_samples <- numeric(0)
    t_samples <- get_nhpp_realization(alpha0,alpha1,alpha2)
    arrivals[i,]<- c(t_samples,rep(NA, col-length(t_samples)))#列数をそろえないと matrix に保存できないため
    c_samples <- hist(arrivals[i,], plot=F, breaks = seq(0, tn, by = 1))$counts
    counts[i,] <- c_samples
  }
  
  #時間データ出力
  write.table(arrivals,paste0(dirname(rstudioapi::getSourceEditorContext()$path),directory_path,"BS-RGS-",estimationMethod,"-",thresholdRule,"-",thresholdValue,"(J=",roopj,").dat"), append=F, quote=F, col.names=F, row.names=F)#出力ファイル名
  #個数データ出力
  write.table(counts,paste0(dirname(rstudioapi::getSourceEditorContext()$path),directory_path,"BSCount-RGS-",estimationMethod,"-",thresholdRule,"-",thresholdValue,"(J=",roopj,").dat"), append=F, quote=F, col.names=F, row.names=F)#出力ファイル名
}