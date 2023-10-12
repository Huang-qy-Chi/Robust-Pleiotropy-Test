#Function of MM-Estimate
rhox = function(x, c = 3.45){   #bisquare function rho
  if(abs(x)<=c){
    rhox1 = 1-(1-(x/c)^(2))^(3)
  }else{
    rhox1 = 1
  }
  return(rhox1)
}

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

#From Maronna(2010):
#q=1,c=1.56;q=2,c=2.66;q=3,c=3.45;q=4,c=4.10;q=5,c=4.65;q=10,c=6.77
#use the following function cacu_c and also obtain a proper c
#Choose the c value
cacu_c = function(q, Sigma = diag(q)){
  c0 = 1:1000/100
  b=0.5
  u = rmvn(1000,mu=rep(0,q),Sigma = Sigma)
  Mr = sqrt(diag(u%*%t(u)))
  for(i in 1:1000){
    c = c0[i]  #use the S-estimator to calculate the scale c
    V1 = (1/1000)*sum(sapply(Mr/c,rhox,c=1))-b
    if(V1<0){
      break
    }
  }
  return(c)
}

#MM-Estimate by Maronna(2010)
#The first function needs to adjust parameters , but when using in simulations, it may be faster.
MM_estimate1X = function(y, x, B0, Sigma0, err = 0.00001){
  #This algorithm mainly base on Kudraszow and Marronna(2010),section 7
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  #if(det(Sigma0)==0){    #avoid singular
  #Sigma0 = Sigma0+err*diag(q)
  #}
  Gamma0 = Sigma0/det(Sigma0)^(1/q) 
  #det(Gamma0)=1
  b = 0.5 #0.5BDP
  ##Step1:M-estimate for scale 
  s0 = (1:1000)/100    #for simplicity we assume s ranges from 0 to 10, which can be adjusted
  #calculate the Mahalanobis distance without scale，using for estimating the scale parameter s
  res = y-x%*%B0
  Mr = sqrt(diag(res%*%solve(Gamma0)%*%t(res)))

  #Step 1: Screening the scale s
  for(i in 1:1000){
    s = s0[i]  #use S-estimator to estimate the scale s
    V1 = (1/n)*sum(sapply(Mr/s,rhox))-b
    if(V1<0){
      break
    }
  }
  
  ##Step2: Estimate of B by WLS
  Sigma0 = s^(2)*Gamma0
  Mr = sqrt(diag(res%*%solve(Sigma0)%*%t(res)))
  weight0 = sapply(Mr,w)
  A2 = t(x)%*%diag(weight0)%*%x
  B2 = t(x)%*%diag(weight0)%*%y
  B1 = solve(A2)%*%B2
  u = y-x%*%B1
  
  ##Step3: Estimate of Gamma and finally Sigma
  C1 = t(u)%*%diag(weight0)%*%u
  Gamma1 = C1/det(C1)^(1/q)
  Sigma1 = s^(2)*Gamma1
  
  ##Step4: Stop until Mahalanobis distance converges
  loop = 1
  Loop = 1000 #max repeat times
  index = 0  #whether the algorithm converge
  while(loop<Loop){
    B0 = B1
    Mr0 = Mr
    #repeat
    res = y-x%*%B0
    #Sigma0 = s^(2)*Gamma0
    Sigma0=Sigma1
    Mr = sqrt(diag(res%*%solve(Sigma0)%*%t(res)))
    
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

#The second function does not need to adjust parameters, using for real data analysis is proper, but it also needs initial values which are affine-equivalent(eg. LTS)
MM_estimate2X = function(y, x, B0, Sigma0, err = 0.00001){
  #This algorithm mainly base on Kudraszow and Marronna(2010),section 7
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  
  cacu_c = function(q, Sigma){ #c ranges from 0 to 10, can be adjusted.
  c0 = 1:1000/100
  b=0.5
  u = rmvn(1000,mu=rep(0,q),Sigma = Sigma)
  Mr = sqrt(diag(u%*%t(u)))
  for(i in 1:1000){
    c = c0[i]  #use S-estimator to calculate the scale s
    V1 = (1/1000)*sum(sapply(Mr/c,rhox,c=1))-b
      if(V1<0){
        break
      }
    }
    return(c)
  }

  c = cacu_c(q, Sigma0)

  rhox = function(x, c){   #bisquare function rho
  if(abs(x)<=c){
    rhox1 = 1-(1-(x/c)^(2))^(3)
    }else{
    rhox1 = 1
    }
  return(rhox1)
  }

  phi = function(x, c){  #first derivative of rho
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

  dw = function(x, c){
    if(abs(x)>=c){
      dw1 = 0
    }else{
      dw1 = -(24/c^(4))*x*(1-(x/c)^(2))
    }
    return(dw1)
  }


  #if(det(Sigma0)==0){    #avoid singular
  #Sigma0 = Sigma0+err*diag(q)
  #}
  Gamma0 = Sigma0/det(Sigma0)^(1/q) 
  #det(Gamma0)=1
  b = 0.5 #0.5BDP
  ##Step1:M-estimate for scale 
  s0 = (1:1000)/100    #for simplicity we assume s ranges from 0 to 10, which can be adjusted
  #calculate the Mahalanobis distance without scale，using for estimating the scale parameter s
  res = y-x%*%B0
  Mr = sqrt(diag(res%*%solve(Gamma0)%*%t(res)))

  
  
  #Step 1: Screening the scale s
  for(i in 1:1000){
    s = s0[i]  #use S-estimator to estimate the scale s
    V1 = (1/n)*sum(sapply(Mr/s,rhox))-b
    if(V1<0){
      break
    }
  }
  
  
  
  ##Step2: Estimate of B by WLS
  Sigma0 = s^(2)*Gamma0
  Mr = sqrt(diag(res%*%solve(Sigma0)%*%t(res)))
  weight0 = sapply(Mr,w)
  A2 = t(x)%*%diag(weight0)%*%x
  B2 = t(x)%*%diag(weight0)%*%y
  B1 = solve(A2)%*%B2
  u = y-x%*%B1
  
  ##Step3: Estimate of Gamma and finally Sigma
  C1 = t(u)%*%diag(weight0)%*%u
  Gamma1 = C1/det(C1)^(1/q)
  Sigma1 = s^(2)*Gamma1
  
  ##Step4: Stop until Mahalanobis distance converges
  loop = 1
  Loop = 1000 #max repeat times
  index = 0  #whether the algorithm converge
  while(loop<Loop){
    #Gamma0 = Gamma1
    #C=solve(Gamma0)
    #u = y-x%*%B0
    B0 = B1
    Mr0 = Mr
    #repeat
    res = y-x%*%B0
    #Sigma0 = s^(2)*Gamma0
    Sigma0=Sigma1
    Mr = sqrt(diag(res%*%solve(Sigma0)%*%t(res)))
    
    
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

