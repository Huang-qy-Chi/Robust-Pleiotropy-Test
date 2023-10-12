#Our proposed robust integrated test: for simulation and may be faster
EHT_PFCX = function(y, x, B0, Sigma0, s, alpha = 0.05){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  u = y-x%*%B0
  invS = solve(Sigma0)
  
  ##Step1:calculate the first and the second derivates
  ###Mahalanobis distance
  Mr = sqrt(diag((u)%*%invS%*%t(u)))
  
  ###the first derivate, only for x is one-dimensional according to the definition of genetic pleiotropy
  #u1 = invS%*%t(u)%*%diag(as.numeric(x))
  u1 = invS%*%t(u)%*%diag(as.numeric(x))
  Phir = u1%*%diag(sapply(Mr,w))
  M = cov(t(Phir))
  
  ###the second derivate and the asymptotic variance of parameter B
  Lambda=(u1)%*%diag(sapply(Mr,dw)/Mr)%*%t(u1)+as.numeric(t(x)%*%diag(sapply(Mr,w))%*%x)*invS
  Lambda1 = Lambda/n
  invLa = solve(Lambda1)
  asyS = invLa%*%M%*%t(invLa)   #B0's asymptotic variance
  
  ##Step2: Test statistics and their values
  BM = B0
  V = diag(q)
  Trecord = n*as.numeric(t(Vec(BM))%*%Lambda1%*%Vec(BM))
  Te = function(m){
    k = m-1
    return(n*as.numeric(t(Vec(BM))%*%t(V[-k,])%*%solve((V[-k,])%*%solve(Lambda1)%*%t(V[-k,]))%*%(V[-k,])%*%Vec(BM)))
  }
  Trecord = c(Trecord, sapply(2:(q+1), Te))
  Tmax = max(Trecord[-1])
  
  ##Step3:calculate the quantiles
  ###calculate the alpha-quantile by Monte-Carlo
  asySq = Msq(asyS,1/2)
  a1 = alpha/(1-alpha)
  quant = rep(0,q+1)
  re0 = matrix(rep(0,1000*(q+1)),nrow = 1000)
  eta = LaplacesDemon::rmvn(1000,rep(0,q),diag(q))
  PV = asySq%*%t(V)%*%solve((V)%*%invLa%*%t(V))%*%(V)%*%asySq
  re0[,1] = diag(eta%*%PV%*%t(eta))  
  quant[1] = (sort(re0[,1])[floor(1000*(1-a1))])
  for(k in 1:q){
    PV = asySq%*%t(V[-k,])%*%solve((V[-k,])%*%invLa%*%t(V[-k,]))%*%(V[-k,])%*%asySq
    re0[,k+1] = diag(eta%*%PV%*%t(eta))
    quant[k+1] = (sort(re0[,k+1])[floor(1000*(1-alpha))])
  } 
  
  ###the value of c1 in our test: a faster way than the paper proposed
  z = which.max(Trecord[-1])
  a2 = a1 - alpha
  c1 = sort(re0[,z])[floor(1000*(1-a2))]

  ##Step4: the decision function: if the decision function equals to 1, we reject H0.
  Te1 = c(Trecord[1]>quant[1],Tmax<c1) 
  Te2 = c(Trecord[-1]>quant[-1])
  index = max(c(sum(Te1)==2,sum(Te2)==q))
  
  return(list(Trecord = Trecord, quant = quant, Test = c(Te1[1],Te2),index = index))
}


#############################################################################################################
#calculate c for our integrated test
comb_c <- function(y, x, B0, Sigma0, alpha=0.05, B=1000){
  q = ncol(y)
  p = ncol(x)
  n = nrow(y)
  B0 = as.matrix(B0)
  x = as.matrix(x)
  Tk= matrix(rep(0,B*(q+1)),ncol = (q+1))
  u = y-x%*%B0[2,]-B0[1,]
  #Sigma0 = Sigma0/det(Sigma0)^(1/q)
  invS = solve(Sigma0)
  ##Step1:calculate the first and the second derivates
  ###Mahalanobis distance
  Mr = sqrt(diag((u)%*%invS%*%t(u)))
  
  ###he first derivate, only for x is one-dimensional according to the definition of genetic pleiotropy
  u1 = invS%*%t(u)%*%diag(as.numeric(x))
  Phir = u1%*%diag(sapply(Mr,w))
  #rowMeans(Phir)
  M = cov(t(Phir))
  
  ###the second derivate and the asymptotic variance of parameter B
  Lambda=(u1)%*%diag(sapply(Mr,dw)/Mr)%*%t(u1)+as.numeric(t(x)%*%diag(sapply(Mr,w))%*%x)*invS
  Lambda1 = Lambda/n
  invLa = solve(Lambda1)
  #eigen(Lambda1)
  asyS = invLa%*%M%*%t(invLa)   #B0's asymptotic variance
  #eigen(asyS)
  
  
  
  ##Step2: Test statistics and their values
  BM = B0
  V = diag(q)
  asySq = Msq(asyS,1/2)
  quant = rep(0,q+1)
  eta = LaplacesDemon::rmvn(B,rep(0,q),diag(q))
  alphas = alpha/(1-alpha)
  PV = asySq%*%t(V)%*%solve((V)%*%invLa%*%t(V))%*%(V)%*%asySq
  Tk[,1] = diag(eta%*%PV%*%t(eta))
  
  for(k in 1:q){
    PV = asySq%*%t(V[-k,])%*%solve((V[-k,])%*%invLa%*%t(V[-k,]))%*%(V[-k,])%*%asySq
    #qu = function(s){
    #return(t(eta[s,])%*%PV%*%eta[s,])
    #}
    #quant[k+1] = (sort(sapply(1:1000,qu))[floor(1000*(1-alpha))])
    Tk[,k+1] = diag(eta%*%PV%*%t(eta))
  } 
  quant[1] = sort(Tk[,1])[floor(B*(1-alphas))]
  for(k in 1:q){
    quant[k+1] = sort(Tk[,k+1])[floor(B*(1-alpha))]
  }
  Tmx = max(Tk[,1])
  div = (Tmx - quant[1])/1000
  grid = seq(quant[1],Tmx,div)
  wq = rep(0,1000)
  for(v in 1:1000){
    c = grid[v]
    index3 = rep(0,B)
    for(d in 1:B){
      Tmax = max(Tk[d,-1])
      Te1 = c(Tk[d,1]>quant[1],Tmax<c) 
      Te2 = c(Tk[d,-1]>quant[-1])
      index3[d] = max(c(sum(Te1)==2,sum(Te2)==q))
    }
    wq[v] = mean(index3)
  }
  asd = which.min(abs(wq-alpha))
  c_res = grid[asd]
  return(c = c_res)
}

#############################################################################################
#Integrated test for real data used, may be slow but clear to use.
EHT_PFC1X = function(y, x, B0, Sigma0, s, alpha = 0.05, B=1000){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  u = y-x%*%B0
  #Sigma0 = Sigma0/det(Sigma0)^(1/q)
  invS = solve(Sigma0)
  
  ##Step1:calculate the first and the second derivates
  ###Mahalanobis distance
  Mr = sqrt(diag((u)%*%invS%*%t(u)))
  
  ###the first derivate, only for x is one-dimensional according to the definition of genetic pleiotropy
  u1 = invS%*%t(u)%*%diag(as.numeric(x))
  Phir = u1%*%diag(sapply(Mr,w))
  M = cov(t(Phir))
  
  ###the second derivate and the asymptotic variance of parameter B
  Lambda=(u1)%*%diag(sapply(Mr,dw)/Mr)%*%t(u1)+as.numeric(t(x)%*%diag(sapply(Mr,w))%*%x)*invS
  Lambda1 = Lambda/n
  invLa = solve(Lambda1)
  #eigen(Lambda1)
  asyS = invLa%*%M%*%t(invLa)#B0渐近方差
  #eigen(asyS)
  
  ##Step2: Test statistics and their values
  BM = B0
  V = diag(q)
  Trecord = n*as.numeric(t(Vec(BM))%*%Lambda1%*%Vec(BM))
  Te = function(m){
    k = m-1
    return(n*as.numeric(t(Vec(BM))%*%t(V[-k,])%*%solve((V[-k,])%*%solve(Lambda1)%*%t(V[-k,]))%*%(V[-k,])%*%Vec(BM)))
  }
  Trecord = c(Trecord, sapply(2:(q+1), Te))
  Tmax = max(Trecord[-1])
  
  ##Step3:calculate the quantiles
  ###calculate the alpha-quantile by Monte-Carlo
  asySq = Msq(asyS,1/2)
  a1 = alpha/(1-alpha)
  quant = rep(0,q+1)
  re0 = matrix(rep(0,1000*(q+1)),nrow = 1000)
  eta = LaplacesDemon::rmvn(1000,rep(0,q),diag(q))
  PV = asySq%*%t(V)%*%solve((V)%*%invLa%*%t(V))%*%(V)%*%asySq
  re0[,1] = diag(eta%*%PV%*%t(eta))
  
  quant[1] = (sort(re0[,1])[floor(1000*(1-a1))])
  for(k in 1:q){
    PV = asySq%*%t(V[-k,])%*%solve((V[-k,])%*%invLa%*%t(V[-k,]))%*%(V[-k,])%*%asySq
    re0[,k+1] = diag(eta%*%PV%*%t(eta))
    quant[k+1] = (sort(re0[,k+1])[floor(1000*(1-alpha))])
  } 
  
  ### the value of c1 in our proposed integrated test
  c = comb_c(y, x, B0, Sigma0, alpha, B)

  
  ##Step4:the decision function: if the decision function equals to 1, we reject H0.
  Te1 = c(Trecord[1]>quant[1],Tmax<c)
  Te2 = c(Trecord[-1]>quant[-1])
  
  #q+1-stage tests are all refused
  index = max(c(sum(Te1)==2,sum(Te2)==q))
  if(index==1){
     print("There exists genetic pleiotropy.")
    }else{
     print("There does not exist genetic pleiotropy.")
    }
  return(list(Trecord = Trecord, quant = quant, Test = c(Te1[1],Te2),index = index))
}
