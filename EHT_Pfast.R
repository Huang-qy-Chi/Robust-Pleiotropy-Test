
#EHT
EHT_PF1 = function(y, x, B0, Sigma0, alpha = 0.05){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  u = y-x%*%B0
  Sigma0 = Sigma0/det(Sigma0)^(1/q)
  invS = solve(Sigma0)
  ##Step1:计算一阶二阶导数，并计算一阶导协方差
  ###马氏距离
  Mr = sqrt(diag((u)%*%invS%*%t(u)))
  
  ###一阶导，仅针对x为一维的情形
  u1 = invS%*%t(u)%*%diag(as.numeric(x))
  Phir = u1%*%diag(sapply(Mr,w))
  #rowMeans(Phir)
  M = cov(t(Phir))
  
  #Phir = w(Mr[1])*Vec(solve(Sigma0)%*%u[1,]%*%(x[1,]),byrow = F)
  #for(h in 2:n){
  #Phir = cbind(Phir, w(Mr[h])*Vec(solve(Sigma0)%*%u[h,]%*%(x[h,]),byrow = F))
  #}
  #Phir-sapply(1:n,ph)
  #M = cov(t(Phir))
  
  ###二阶导：极小值存在，理论上二阶导正定
  Lambda=(u1)%*%diag(sapply(Mr,dw)/Mr)%*%t(u1)+as.numeric(t(x)%*%diag(sapply(Mr,w))%*%x)*invS
  Lambda1 = Lambda/n
  invLa = solve(Lambda1)
  #eigen(Lambda1)
  asyS = invLa%*%M%*%t(invLa)#B0渐近方差
  #eigen(asyS)
  
  
  
  ##Step2:构成统计量并计算统计量的值
  BM = B0
  V = diag(q)
  Trecord = n*as.numeric(t(Vec(BM))%*%Lambda1%*%Vec(BM))
  Te = function(m){
    k = m-1
    return(n*as.numeric(t(Vec(BM))%*%t(V[-k,])%*%solve((V[-k,])%*%solve(Lambda1)%*%t(V[-k,]))%*%(V[-k,])%*%Vec(BM)))
  }
  Trecord = c(Trecord, sapply(2:(q+1), Te))
  
  ##Step3:计算分位数
  ###批量求解特征值
  ###求解alpha分位数：避免重复计算矩阵
  asySq = Msq(asyS,1/2)
  quant = rep(0,q+1)
  eta = LaplacesDemon::rmvn(1000,rep(0,q),diag(q))
  PV = asySq%*%t(V)%*%solve((V)%*%invLa%*%t(V))%*%(V)%*%asySq
  #qu = function(k){
    #return(t(eta[k,])%*%PV%*%(eta[k,]))
  #}
  #quant[1] = (sort(sapply(1:1000,qu))[floor(1000*(1-alpha))])
  quant[1] = sort(diag(eta%*%PV%*%t(eta)))[floor(1000*(1-alpha))]
  #k=2
  for(k in 1:q){
    PV = asySq%*%t(V[-k,])%*%solve((V[-k,])%*%invLa%*%t(V[-k,]))%*%(V[-k,])%*%asySq
    #qu = function(s){
      #return(t(eta[s,])%*%PV%*%eta[s,])
    #}
    #quant[k+1] = (sort(sapply(1:1000,qu))[floor(1000*(1-alpha))])
    quant[k+1] = sort(diag(eta%*%PV%*%t(eta)))[floor(1000*(1-alpha))]
  } 
  
  #
  
  
  ##Step4:找出显著不为0的分量
  index = 0
  Te = c(Trecord>quant) #改用逻辑判断代替循环
  if(sum(Te)==(q+1)){   #q+1-stage tests are all refused
    index = 1
  }
  return(list(Trecord = Trecord, quant = quant, Test = Te,index = index))
}





EHT_PFC = function(y, x, B0, Sigma0, alpha = 0.05){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  u = y-x%*%B0
  Sigma0 = Sigma0/det(Sigma0)^(1/q)
  invS = solve(Sigma0)
  ##Step1:计算一阶二阶导数，并计算一阶导协方差
  ###马氏距离
  Mr = sqrt(diag((u)%*%invS%*%t(u)))
  
  ###一阶导，仅针对x为一维的情形
  #u1 = invS%*%t(u)%*%diag(as.numeric(x))
  u1 = invS%*%t(u)%*%diag(as.numeric(x))
  Phir = u1%*%diag(sapply(Mr,w))
  #mean(Phir)
  M = cov(t(Phir))
  
  #Phir = w(Mr[1])*Vec(solve(Sigma0)%*%u[1,]%*%(x[1,]),byrow = F)
  #for(h in 2:n){
  #Phir = cbind(Phir, w(Mr[h])*Vec(solve(Sigma0)%*%u[h,]%*%(x[h,]),byrow = F))
  #}
  #Phir-sapply(1:n,ph)
  #M = cov(t(Phir))
  
  ###二阶导：极小值存在，理论上二阶导正定
  Lambda=(u1)%*%diag(sapply(Mr,dw)/Mr)%*%t(u1)+as.numeric(t(x)%*%diag(sapply(Mr,w))%*%x)*invS
  Lambda1 = Lambda/n
  invLa = solve(Lambda1)
  #eigen(Lambda1)
  asyS = invLa%*%M%*%t(invLa)#B0渐近方差
  #eigen(asyS)
  
  
  
  ##Step2:构成统计量并计算统计量的值
  BM = B0
  V = diag(q)
  Trecord = n*as.numeric(t(Vec(BM))%*%Lambda1%*%Vec(BM))
  Te = function(m){
    k = m-1
    return(n*as.numeric(t(Vec(BM))%*%t(V[-k,])%*%solve((V[-k,])%*%solve(Lambda1)%*%t(V[-k,]))%*%(V[-k,])%*%Vec(BM)))
  }
  Trecord = c(Trecord, sapply(2:(q+1), Te))
  Tmax = max(Trecord[-1])
  ##Step3:计算分位数
  ###批量求解特征值
  ###求解alpha分位数：避免重复计算矩阵
  asySq = Msq(asyS,1/2)
  a1 = alpha/(1-alpha)
  quant = rep(0,q+1)
  re0 = matrix(rep(0,1000*(q+1)),nrow = 1000)
  eta = LaplacesDemon::rmvn(1000,rep(0,q),diag(q))
  PV = asySq%*%t(V)%*%solve((V)%*%invLa%*%t(V))%*%(V)%*%asySq
  re0[,1] = diag(eta%*%PV%*%t(eta))
  
  quant[1] = (sort(re0[,1])[floor(1000*(1-a1))])
  #k=2
  for(k in 1:q){
    PV = asySq%*%t(V[-k,])%*%solve((V[-k,])%*%invLa%*%t(V[-k,]))%*%(V[-k,])%*%asySq
    re0[,k+1] = diag(eta%*%PV%*%t(eta))
    quant[k+1] = (sort(re0[,k+1])[floor(1000*(1-alpha))])
  } 
  
  #求解c的经验值
  z = which.max(Trecord[-1])
  #a2 = alpha^(2)/(1-alpha)
  a2 = a1 - alpha
  c = sort(re0[,z])[floor(1000*(1-a2))]
  #plot(density(sort(re0[,z])))
    
  
  
  
  ##Step4:找出显著不为0的分量
  Te1 = c(Trecord[1]>quant[1],Tmax<c) #改用逻辑判断代替循环
  Te2 = c(Trecord[-1]>quant[-1])
  
     #q+1-stage tests are all refused
  index = max(c(sum(Te1)==2,sum(Te2)==q))
  
  return(list(Trecord = Trecord, quant = quant, Test = c(Te1[1],Te2),index = index))
}



#带截距项
EHT_PFCM = function(y, x, B0, Sigma0, alpha = 0.05){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  B0 = as.matrix(B0)
  x = as.matrix(x)
  u = y-x%*%B0[2,]-B0[1,]
  Sigma0 = Sigma0/det(Sigma0)^(1/q)
  invS = solve(Sigma0)
  ##Step1:计算一阶二阶导数，并计算一阶导协方差
  ###马氏距离
  Mr = sqrt(diag((u)%*%invS%*%t(u)))
  
  x0 = x-mean(x)
  y0 = y - matrix(rep((mean(x)%*%B0[2,]),n),nrow = n)
  x = x0
  y = y0
  ###一阶导，仅针对x为一维的情形
  #u1 = invS%*%t(u)%*%diag(as.numeric(x))
  u1 = invS%*%t(u)%*%diag(as.numeric(x))
  Phir = u1%*%diag(sapply(Mr,w))
  #mean(Phir)
  M = cov(t(Phir))
  
  #Phir = w(Mr[1])*Vec(solve(Sigma0)%*%u[1,]%*%(x[1,]),byrow = F)
  #for(h in 2:n){
  #Phir = cbind(Phir, w(Mr[h])*Vec(solve(Sigma0)%*%u[h,]%*%(x[h,]),byrow = F))
  #}
  #Phir-sapply(1:n,ph)
  #M = cov(t(Phir))
  
  ###二阶导：极小值存在，理论上二阶导正定
  Lambda=(u1)%*%diag(sapply(Mr,dw)/Mr)%*%t(u1)+as.numeric(t(x)%*%diag(sapply(Mr,w))%*%x)*invS
  Lambda1 = Lambda/n
  invLa = solve(Lambda1)
  #eigen(Lambda1)
  asyS = invLa%*%M%*%t(invLa)#B0渐近方差
  #eigen(asyS)
  
  
  
  ##Step2:构成统计量并计算统计量的值
  BM = B0
  V = diag(q)
  Trecord = n*as.numeric(t((BM[2,]))%*%Lambda1%*%(BM[2,]))
  Te = function(m){
    k = m-1
    return(n*as.numeric(t((BM[2,]))%*%t(V[-k,])%*%solve((V[-k,])%*%solve(Lambda1)%*%t(V[-k,]))%*%(V[-k,])%*%(BM[2,])))
  }
  Trecord = c(Trecord, sapply(2:(q+1), Te))
  Tmax = max(Trecord[-1])
  ##Step3:计算分位数
  ###批量求解特征值
  ###求解alpha分位数：避免重复计算矩阵
  asySq = Msq(asyS,1/2)
  a1 = alpha/(1-alpha)
  quant = rep(0,q+1)
  re0 = matrix(rep(0,1000*(q+1)),nrow = 1000)
  eta = LaplacesDemon::rmvn(1000,rep(0,q),diag(q))
  PV = asySq%*%t(V)%*%solve((V)%*%invLa%*%t(V))%*%(V)%*%asySq
  re0[,1] = diag(eta%*%PV%*%t(eta))
  
  quant[1] = (sort(re0[,1])[floor(1000*(1-a1))])
  #k=2
  for(k in 1:q){
    PV = asySq%*%t(V[-k,])%*%solve((V[-k,])%*%invLa%*%t(V[-k,]))%*%(V[-k,])%*%asySq
    re0[,k+1] = diag(eta%*%PV%*%t(eta))
    quant[k+1] = (sort(re0[,k+1])[floor(1000*(1-alpha))])
  } 
  
  #求解c的经验值
  z = which.max(Trecord[-1])
  #a2 = alpha^(2)/(1-alpha)
  a2 = a1 - alpha
  c = sort(re0[,z])[floor(1000*(1-a2))]
  #plot(density(sort(re0[,z])))
  
  
  
  
  ##Step4:找出显著不为0的分量
  Te1 = c(Trecord[1]>quant[1],Tmax<c) #改用逻辑判断代替循环
  Te2 = c(Trecord[-1]>quant[-1])
  
  #q+1-stage tests are all refused
  index = max(c(sum(Te1)==2,sum(Te2)==q))
  
  return(list(Trecord = Trecord, quant = quant, Test = c(Te1[1],Te2),index = index))
}

#带截距项，q+1检验
EHT_PF1M = function(y, x, B0, Sigma0, alpha = 0.05){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  B0 = as.matrix(B0)
  x = as.matrix(x)
  u = y-x%*%B0[2,]-B0[1,]
  Sigma0 = Sigma0/det(Sigma0)^(1/q)
  invS = solve(Sigma0)
  ##Step1:计算一阶二阶导数，并计算一阶导协方差
  ###马氏距离
  Mr = sqrt(diag((u)%*%invS%*%t(u)))
  
  x0 = x-mean(x)
  y0 = y - matrix(rep((mean(x)%*%B0[2,]),n),nrow = n)
  x = x0
  y = y0
  ###一阶导，仅针对x为一维的情形
  #u1 = invS%*%t(u)%*%diag(as.numeric(x))
  u1 = invS%*%t(u)%*%diag(as.numeric(x))
  Phir = u1%*%diag(sapply(Mr,w))
  #mean(Phir)
  M = cov(t(Phir))
  
  #Phir = w(Mr[1])*Vec(solve(Sigma0)%*%u[1,]%*%(x[1,]),byrow = F)
  #for(h in 2:n){
  #Phir = cbind(Phir, w(Mr[h])*Vec(solve(Sigma0)%*%u[h,]%*%(x[h,]),byrow = F))
  #}
  #Phir-sapply(1:n,ph)
  #M = cov(t(Phir))
  
  ###二阶导：极小值存在，理论上二阶导正定
  Lambda=(u1)%*%diag(sapply(Mr,dw)/Mr)%*%t(u1)+as.numeric(t(x)%*%diag(sapply(Mr,w))%*%x)*invS
  Lambda1 = Lambda/n
  invLa = solve(Lambda1)
  #eigen(Lambda1)
  asyS = invLa%*%M%*%t(invLa)#B0渐近方差
  #eigen(asyS)
  
  
  
  ##Step2:构成统计量并计算统计量的值
  BM = B0
  V = diag(q)
  Trecord = n*as.numeric(t((BM[2,]))%*%Lambda1%*%(BM[2,]))
  Te = function(m){
    k = m-1
    return(n*as.numeric(t((BM[2,]))%*%t(V[-k,])%*%solve((V[-k,])%*%solve(Lambda1)%*%t(V[-k,]))%*%(V[-k,])%*%(BM[2,])))
  }
  Trecord = c(Trecord, sapply(2:(q+1), Te))
  Tmax = max(Trecord[-1])
  ##Step3:计算分位数
  ###批量求解特征值
  ###求解alpha分位数：避免重复计算矩阵
  asySq = Msq(asyS,1/2)
  a1 = alpha/(1-alpha)
  quant = rep(0,q+1)
  re0 = matrix(rep(0,1000*(q+1)),nrow = 1000)
  eta = LaplacesDemon::rmvn(1000,rep(0,q),diag(q))
  PV = asySq%*%t(V)%*%solve((V)%*%invLa%*%t(V))%*%(V)%*%asySq
  re0[,1] = diag(eta%*%PV%*%t(eta))
  
  quant[1] = (sort(re0[,1])[floor(1000*(1-a1))])
  #k=2
  for(k in 1:q){
    PV = asySq%*%t(V[-k,])%*%solve((V[-k,])%*%invLa%*%t(V[-k,]))%*%(V[-k,])%*%asySq
    re0[,k+1] = diag(eta%*%PV%*%t(eta))
    quant[k+1] = (sort(re0[,k+1])[floor(1000*(1-alpha))])
  } 

  ##Step4:找出显著不为0的分量
  Te1 = c(Trecord>quant) #改用逻辑判断代替循环
  #q+1-stage tests are all refused
  index = sum(Te1)==q+1
  return(list(Trecord = Trecord, quant = quant, Test = Te1,index = index))
}

