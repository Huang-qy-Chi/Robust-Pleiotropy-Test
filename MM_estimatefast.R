#Function of MM-Estimate
rhox = function(x, c = 3.45){   #bisquare function rho
  if(abs(x)<=c){
    rhox1 = 1-(1-(x/c)^(2))^(3)
  }else{
    rhox1 = 1
  }
  return(rhox1)
}
#From Maronna(2010):
#q=1,c=1.56;q=2,c=2.66;q=3,c=3.45;q=4,c=4.10;q=5,c=4.65;q=10,c=6.77

phi = function(x,c = 3.45){  #first derivative of rho
  if(abs(x)>=c){
    phix = 0
  }else{
    phix = (6/c^(2))*x*(1-(x/c)^(2))^(2)
  }
  return(phix)
}

w = function(u){  #weight function
  return(phi(u)/u)
}

dw = function(x, c = 3.45){
  if(abs(x)>=c){
    dw1 = 0
  }else{
    dw1 = -(24/c^(4))*x*(1-(x/c)^(2))
  }
  return(dw1)
}

#Choose the c value
cacu_c = function(q){
  c0 = 1:1000/100
  b=0.5
  u = rmvn(1000,mu=rep(0,q),Sigma = diag(q))
  Mr = sqrt(diag(u%*%t(u)))
  for(i in 1:1000){
    c = c0[i]  #使用S估计计算scale
    V1 = (1/1000)*sum(sapply(Mr/c,rhox,c=1))-b
    if(V1<0){
      break
    }
  }
  return(c)
}
#fg = rep(0,1000)
#for(i in 1:1000){
  #fg[i] = cacu_c(6)
#}
#mean(fg)

#MM-Estimate by Maronna(2010)
MM_estimate1 = function(y, x, B0, Sigma0, err = 0.00001){
  #This algorithm mainly base on Kudraszow and Marronna(2010),section 7
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  #if(det(Sigma0)==0){
    #Sigma0 = Sigma0+err*diag(q)
  #}
  Gamma0 = Sigma0/det(Sigma0)^(1/q) #det(Gamma0)=1
  b = 0.5 #0.5BDP
  ##Step1:M-estimate for scale 
  s0 = (1:1000)/100 #初值估计精度
  #计算去除scale的马氏距离，用于估计scale与beta
  res = y-x%*%B0
  Mr = sqrt(diag(res%*%solve(Gamma0)%*%t(res)))
  #Mr = as.numeric(sqrt((y[1,]-x[1,]%*%B0)%*%solve(Gamma0)%*%t(y[1,]-x[1,]%*%B0)))
  #for(j in 2:n){
  #Mr = c(Mr, as.numeric(sqrt((y[j,]-x[j,]%*%B0)%*%solve(Gamma0)%*%t(y[j,]-x[j,]%*%B0))))
  #}
  #Mr-Mr1summary(Mr)
  #system.time({
    #for(i in 1:10000){
      #seq(1:1000)/100
      #(1:1000)/100
    #}
  #})
  
  #栅格法近似计算scale
  for(i in 1:1000){
    s = s0[i]  #使用S估计计算scale
    V1 = (1/n)*sum(sapply(Mr/s,rhox))-b
    if(V1<0){
      break
    }
  }
  
  
  
  ##Step2: Estimate of B by WLS
  weight0 = sapply(Mr,w)
  A2 = t(x)%*%diag(weight0)%*%x
  B2 = t(x)%*%diag(weight0)%*%y
  B1 = solve(A2)%*%B2
  u = y-x%*%B1
  
  ##Step3: Estimate of Gamma
  C1 = t(u)%*%diag(weight0)%*%u
  Gamma1 = C1/det(C1)^(1/q)
  Sigma1 = s^(2)*Gamma1
  
  ##Step4: Convergence
  loop = 1
  Loop = 1000
  index = 0
  while(loop<Loop){
    Gamma0 = Gamma1
    #C=solve(Gamma0)
    #u = y-x%*%B0
    B0 = B1
    Mr0 = Mr
    #重复迭代
    res = y-x%*%B0
    Mr = sqrt(diag(res%*%solve(Gamma0)%*%t(res)))
    #Mr = as.numeric(sqrt((y[1,]-x[1,]%*%B0)%*%solve(Gamma0)%*%t(y[1,]-x[1,]%*%B0)))
    #for(j in 2:n){
    #Mr = c(Mr, as.numeric(sqrt((y[j,]-x[j,]%*%B0)%*%solve(Gamma0)%*%t(y[j,]-x[j,]%*%B0))))
    #}
    #Mr-Mr1
    
    ##Step2: Estimate of B by WLS
    weight0 = sapply(Mr,w)
    A2 = t(x)%*%diag(weight0)%*%x
    B2 = t(x)%*%diag(weight0)%*%y
    B1 = solve(A2)%*%B2
    u = y-x%*%B1
    
    ##Step3: Estimate of Gamma
    C1 = t(u)%*%diag(weight0)%*%u
    Gamma1 = C1/det(C1)^(1/q)
    Sigma1 = s^(2)*Gamma1
    if(t(Mr-Mr0)%*%(Mr-Mr0)<err){
      index = 1
      break
    }else{
      #cat("Loop",loop,"times.\r")
      loop = loop + 1
    }
  }
  return(list(Beta = B1, Sigma = Sigma1,s=s, index = index))
}
#MM_estimate1(y,x,B0,Sigma0)


#带截距项
MM_estimate2 = function(y, x, B0, Sigma0, err = 0.00001){
  #This algorithm mainly base on Kudraszow and Marronna(2010),section 7
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  x = cbind(rep(1,n),x)
  #if(det(Sigma0)==0){
  #Sigma0 = Sigma0+err*diag(q)
  #}
  Gamma0 = Sigma0/det(Sigma0)^(1/q) #det(Gamma0)=1
  b = 0.5 #0.5BDP
  ##Step1:M-estimate for scale 
  s0 = (1:1000)/100 #初值估计精度
  #计算去除scale的马氏距离，用于估计scale与beta
  res = y-x%*%B0
  Mr = sqrt(diag(res%*%solve(Gamma0)%*%t(res)))
  #Mr = as.numeric(sqrt((y[1,]-x[1,]%*%B0)%*%solve(Gamma0)%*%t(y[1,]-x[1,]%*%B0)))
  #for(j in 2:n){
  #Mr = c(Mr, as.numeric(sqrt((y[j,]-x[j,]%*%B0)%*%solve(Gamma0)%*%t(y[j,]-x[j,]%*%B0))))
  #}
  #Mr-Mr1summary(Mr)
  #system.time({
  #for(i in 1:10000){
  #seq(1:1000)/100
  #(1:1000)/100
  #}
  #})
  
  #栅格法近似计算scale
  for(i in 1:1000){
    s = s0[i]  #使用S估计计算scale
    V1 = (1/n)*sum(sapply(Mr/s,rhox))-b
    if(V1<0){
      break
    }
  }
  
  
  
  ##Step2: Estimate of B by WLS
  weight0 = sapply(Mr,w)
  A2 = t(x)%*%diag(weight0)%*%x
  B2 = t(x)%*%diag(weight0)%*%y
  B1 = solve(A2)%*%B2
  u = y-x%*%B1
  
  ##Step3: Estimate of Gamma
  C1 = t(u)%*%diag(weight0)%*%u
  Gamma1 = C1/det(C1)^(1/q)
  Sigma1 = s^(2)*Gamma1
  
  ##Step4: Convergence
  loop = 1
  Loop = 1000
  index = 0
  while(loop<Loop){
    Gamma0 = Gamma1
    #C=solve(Gamma0)
    #u = y-x%*%B0
    B0 = B1
    Mr0 = Mr
    #重复迭代
    res = y-x%*%B0
    Mr = sqrt(diag(res%*%solve(Gamma0)%*%t(res)))
    #Mr = as.numeric(sqrt((y[1,]-x[1,]%*%B0)%*%solve(Gamma0)%*%t(y[1,]-x[1,]%*%B0)))
    #for(j in 2:n){
    #Mr = c(Mr, as.numeric(sqrt((y[j,]-x[j,]%*%B0)%*%solve(Gamma0)%*%t(y[j,]-x[j,]%*%B0))))
    #}
    #Mr-Mr1
    
    ##Step2: Estimate of B by WLS
    weight0 = sapply(Mr,w)
    A2 = t(x)%*%diag(weight0)%*%x
    B2 = t(x)%*%diag(weight0)%*%y
    B1 = solve(A2)%*%B2
    u = y-x%*%B1
    
    ##Step3: Estimate of Gamma
    C1 = t(u)%*%diag(weight0)%*%u
    Gamma1 = C1/det(C1)^(1/q)
    Sigma1 = s^(2)*Gamma1
    if(t(Mr-Mr0)%*%(Mr-Mr0)<err){
      index = 1
      break
    }else{
      #cat("Loop",loop,"times.\r")
      loop = loop + 1
    }
  }
  return(list(Beta = B1, Sigma = Sigma1,s=s, index = index))
}




