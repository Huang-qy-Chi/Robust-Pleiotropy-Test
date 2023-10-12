#Schaid et al. (2016), Two Stage Test
EHT_MLH = function(y,x,alpha=0.05){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  V = diag(q)
  x = x - mean(x)
  sX = solve(t(x)%*%(x))
  beta = sX%*%(t(x)%*%y)
  Sigma = (1/n)*t(y)%*%(diag(n)-x%*%sX%*%t(x))%*%y
  Omegasq = kronecker(Msq(Sigma,-0.5),diag(n))
  X = kronecker(diag(q),x)
  Y = Vec(y,byrow = F)
  beta1 = solve(t(X)%*%Omegasq%*%(X))%*%(t(X)%*%Omegasq%*%(Y))
  #standalized
  X1 = Omegasq%*%X
  Y1 = Omegasq%*%Y
  V = diag(q)
  beta2 = solve(t(X1)%*%(X1))%*%(t(X1)%*%(Y1))
  HX = t(X1)%*%(X1)
  Trecord = t(beta2)%*%HX%*%beta2
  for(k in 1:q){
    Trecord = c(Trecord,t(beta2)%*%t(V[-k,])%*%solve(V[-k,]%*%solve(HX)%*%t(V[-k,]))%*%(V[-k,])%*%beta2)
  }
  Tmin = min(Trecord[-1])
  #index = 0
  #if(Trecord[1]>qchisq(df = q,1-alpha)){
    #if(Tmin>qchisq(df = q-1,1-alpha)){
      #index = 1
    #}
  #}
  index = (Trecord[1]>qchisq(df = q,1-alpha)&Tmin>qchisq(df = q-1,1-alpha))
  return(list(index = index,beta = beta2, Trecord = Trecord, quantile = c(qchisq(df = q,1-alpha),qchisq(df = q-1,1-alpha))))
}

#Two Stage Test with combined decision function
EHT_MLC = function(y,x,alpha=0.05){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  V = diag(q)
  x = x - mean(x)
  sX = solve(t(x)%*%(x))
  beta = sX%*%(t(x)%*%y)
  Sigma = (1/n)*t(y)%*%(diag(n)-x%*%sX%*%t(x))%*%y
  Omegasq = kronecker(Msq(Sigma,-0.5),diag(n))
  X = kronecker(diag(q),x)
  Y = Vec(y,byrow = F)
  beta1 = solve(t(X)%*%Omegasq%*%(X))%*%(t(X)%*%Omegasq%*%(Y))
  #standalized
  X1 = Omegasq%*%X
  Y1 = Omegasq%*%Y
  V = diag(q)
  beta2 = solve(t(X1)%*%(X1))%*%(t(X1)%*%(Y1))
  HX = t(X1)%*%(X1)
  Trecord = t(beta2)%*%HX%*%beta2
  for(k in 1:q){
    Trecord = c(Trecord,t(beta2)%*%t(V[-k,])%*%solve(V[-k,]%*%solve(HX)%*%t(V[-k,]))%*%(V[-k,])%*%beta2)
  }
  Tmin = min(Trecord[-1])
  Tmax = max(Trecord[-1])
  index = 0
  alpha1 = alpha/(1-alpha)
  alpha2 = alpha1 - alpha
  TE1 = (Trecord[1]>qchisq(df=q,1-alpha1)&Tmax<qchisq(df=q-1,1-alpha2))
  TE2 = (Tmin>qchisq(df=q-1,1-alpha))
  index = max(TE1,TE2)
  return(list(index = index, Trecord = Trecord, quantile = c(qchisq(df = q,1-alpha),qchisq(df = q-1,1-alpha))))
}


#Methods of Wang et al. (2021)
EHT_MLW = function(y,x,alpha=0.05){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  V = diag(q)
  x = x - mean(x)
  sX = solve(t(x)%*%(x))
  beta = sX%*%(t(x)%*%y)
  Sigma = (1/n)*t(y)%*%(diag(n)-x%*%sX%*%t(x))%*%y
  Omegasq = kronecker(Msq(Sigma,-0.5),diag(n))
  X = kronecker(diag(q),x)
  Y = Vec(y,byrow = F)
  #beta1 = solve(t(X)%*%Omegasq%*%(X))%*%(t(X)%*%Omegasq%*%(Y))
  #standalized
  X1 = Omegasq%*%X
  Y1 = Omegasq%*%Y
  V = diag(q)
  beta2 = solve(t(X1)%*%(X1))%*%(t(X1)%*%(Y1))
  HX = t(X1)%*%(X1)
  Trecord = t(beta2)%*%HX%*%beta2
  for(k in 1:q){
    Trecord = c(Trecord,t(beta2)%*%t(V[-k,])%*%solve(V[-k,]%*%solve(HX)%*%t(V[-k,]))%*%(V[-k,])%*%beta2)
  }
  Tmin = min(Trecord[-1])
  #index = 0
  #if(Trecord[1]>qchisq(df = q,1-alpha)){
  #if(Tmin>qchisq(df = q-1,1-alpha)){
  #index = 1
  #}
  #}
  beta_V = rep(0,q)
  for(j in 1:q){
    betastar = slove(HX)%*%T(V[-k,])%*%t(V[-k,])%*%solve(V[-k,]%*%solve(HX)%*%t(V[-k,]))%*%v[-k,]%8%beta2
    beta_V[j] = beta2 - betastar
  }
  
  likh = function(Y1,X1,beta){ #likelihood
    y = Y1
    x = X1
    return(exp(-t(y-x%*%beta)%*%(y-x%*%beta))) #let the weights greater than 0
  }
  Lik = sapply(c(beta2,beta_V),likh,Y1 = Y1,X1 = X1) #likelihoods
  ww = c(1,sqrt(n)*rep(1,q))
  bic = Lik/ww
  sbic = sum(bic)
  w = bic/sbic
  i_max = which.max(w) #the subhypothesis with the hightest probability
  #calculate p-value
  p_value = rep(0,q+1)
  p_value[1]=1-pchisq(Trecord[1],df = q)
  for(k in 1:q){
    p_value[k+1]=1-pchisq(Trecord[1],df = q-1)
  }
  pm = p_value[i_max]
  index = (pm<alpha)
  return(list(index = index,p_value = pm, Trecord = Trecord, quantile = c(qchisq(df = q,1-alpha),qchisq(df = q-1,1-alpha))))
}
