#By Multipulier MCD: using asymptotic chi square
M_LTS2 = function(y, x, err = 0.00001,alpha = 0.05){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)  #适应多元情形
  value = qchisq(alpha, df = q, lower.tail = F )
  beta0 = rep(0.2, q)
  u = y-x%*%t(beta0)
  #计算协方差
  Sigma0 = (1/n)*t(y-x%*%t(beta0))%*%(y-x%*%t(beta0))
  #Sigma0 = diag(q)
  #计算马氏距离
  
  rsq = diag((u)%*%solve(Sigma0)%*%t(u))
  
  w0 = c(rsq<value)  #逻辑值可以直接返回为数
  
  
  #加权最小二乘估计
  A1 = t(x)%*%diag(w0)%*%x
  B1 = t(x)%*%diag(w0)%*%y
  beta1 = solve(A1)%*%B1
  
  #求解协方差Sigma
  res = y - x%*%(beta1)
  Sigma = t(res)%*%diag(w0)%*%res/sum(w0)
  
  loop = 1
  Loop = 1000
  index = 0
  while(loop<Loop){
    beta0 = beta1
    
    rsq = diag((u)%*%solve(Sigma0)%*%t(u))
    
    w0 = c(rsq<value)  #逻辑值可以直接返回为数
    
    
    #加权最小二乘估计
    A1 = t(x)%*%diag(w0)%*%x
    B1 = t(x)%*%diag(w0)%*%y
    beta1 = solve(A1)%*%B1
    
    #求解协方差Sigma
    res = y - x%*%(beta1)
    Sigma = t(res)%*%diag(w0)%*%res/sum(w0)
    if((beta1-beta0)%*%t(beta1-beta0)<err){
      index = 1
      break
    }else{
      #cat("Loop:",loop,"times.\r")
      loop = loop + 1
    }
  }
  return(list(beta = beta1, Sigma = Sigma, index = index))
}

#library(LaplacesDemon)
#system.time({
 # for(i in 1:1000){
 # yh = Gen_data(100,beta = c(1,0,0),Sigma0 = diag(3))
 # y = yh$Y
 # x = yh$X
  # df1 = M_LTS2(y,x)
   #B0 = df1$beta
   #Sigma0 = df1$Sigma
  # MM_estimate1(y,x,B0,Sigma0)
  # cat("Loop:",i,"times.\r")
# }
#})


#带截距项
M_LTS3 = function(y, x, err = 0.00001,alpha = 0.05){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)  #适应多元情形
  x = cbind(rep(1,n),x)
  value = qchisq(alpha, df = q, lower.tail = F )
  beta0 = matrix(rep(0.2, 2*q),ncol=2)
  u = y-x%*%t(beta0)
  #计算协方差
  Sigma0 = (1/n)*t(y-x%*%t(beta0))%*%(y-x%*%t(beta0))
  #Sigma0 = diag(q)
  #计算马氏距离
  
  rsq = diag((u)%*%solve(Sigma0)%*%t(u))
  
  w0 = c(rsq<value)  #逻辑值可以直接返回为数
  
  
  #加权最小二乘估计
  A1 = t(x)%*%diag(w0)%*%x
  B1 = t(x)%*%diag(w0)%*%y
  beta1 = solve(A1)%*%B1
  
  #求解协方差Sigma
  res = y - x%*%(beta1)
  Sigma = t(res)%*%diag(w0)%*%res/sum(w0)
  
  loop = 1
  Loop = 1000
  index = 0
  while(loop<Loop){
    beta0 = beta1
    
    rsq = diag((u)%*%solve(Sigma0)%*%t(u))
    
    w0 = c(rsq<value)  #逻辑值可以直接返回为数
    
    
    #加权最小二乘估计
    A1 = t(x)%*%diag(w0)%*%x
    B1 = t(x)%*%diag(w0)%*%y
    beta1 = solve(A1)%*%B1
    
    #求解协方差Sigma
    res = y - x%*%(beta1)
    Sigma = t(res)%*%diag(w0)%*%res/sum(w0)
    if(sum(diag((beta1-beta0)%*%t(beta1-beta0)))<err){
      index = 1
      break
    }else{
      #cat("Loop:",loop,"times.\r")
      loop = loop + 1
    }
  }
  return(list(beta = beta1, Sigma = Sigma, index = index))
}

