# Function of Data Genetate
library(MASS)
library(LaplacesDemon)


#normal
library(LaplacesDemon)
Gen_data = function(n = 100, beta, Sigma0){
  q = length(beta)
  X = rnorm(n,0,1)
  epsilon = rmvn(n, mu = rep(0,q), Sigma = Sigma0)
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}

#cauchy
Gen_data2 = function(n = 100, beta, S, COR){
  library(LaplacesDemon)
  q = length(beta)
  X = rnorm(n)
  epsiolon = rmvc(n, mu = rep(0,q), S = S)%*%Msq(COR,0.5)
  Y = X%*%t(beta) + epsiolon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}

#mixnormal
Gen_data3 = function(n = 100, beta, mu, Sigma1, Sigma2,ratio = 0.1){
  library(LaplacesDemon)
  q = as.numeric(ncol(Sigma1))
  X = rnorm(n)
  #epsilon = rmvn(n, mu = rep(0,q), Sigma = Sigma1)+rbind(matrix(rep(0,n*(1-ratio)*q),ncol=q),rmvn(n*ratio,mu,Sigma2))
  #epsilon = 0.8*rmvn(n, mu = rep(0,q), Sigma = Sigma1)+0.2*rmvn(n, mu = mu, Sigma = Sigma2)
  w = runif(n)
  epsilon = matrix(rep(0,n*q),ncol = q)
  for(j in 1:n){
    if(w[j]<1-ratio){
      epsilon[j,] = rmvn(1, mu = rep(0,q), Sigma = Sigma1)
    }else{
      epsilon[j,] = rmvn(1, mu = mu, Sigma = Sigma2)
    }
  }
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}

#Question heteroscedasticity
Gen_data4 = function(n = 100, beta, Sigma0){
  library(LaplacesDemon)
  q = as.numeric(ncol(Sigma0))
  X = rnorm(n)
  epsilon = rmvn(1, mu = rep(0,q), Sigma = (1+(X[1]^(2)))*Sigma0)
  for(i in 2:n){
    epsilon = rbind(epsilon, rmvn(1, mu = rep(0,q), Sigma = (1+(X[i]^(2)))*Sigma0))
  }
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}



library(LaplacesDemon)
# real heteroscedasticity
Gen_data5 = function(n = 100, beta, Sigma0){
  q = as.numeric(ncol(Sigma0))
  X = rnorm(n)
  epsilon = (diag(n)+1*diag(X))%*%rmvn(n, mu = rep(0,q), Sigma = Sigma0)
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}



#test for distribution of x
Gen_data6 = function(n = 100, beta, Sigma0){
  q = length(beta)
  #X = rnorm(n,5,1)
  #X = rcauchy(n)
  X = rnorm(n)
  epsilon = rmvn(n, mu = rep(0,q), Sigma = Sigma0)
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}

#t distribution
Gen_data7 = function(n = 100, beta, Sigma0, df = Inf){
  q = length(beta)
  #X = rnorm(n,5,1)
  #X = rcauchy(n)
  X = rnorm(n)
  epsilon = rmvt(n, mu = rep(0,q), S = Sigma0, df)
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}

#Multivarate laplace
Gen_data8 = function(n = 100, beta, Sigma0){
  q = length(beta)
  #X = rnorm(n,5,1)
  #X = rcauchy(n)
  X = rnorm(n)
  epsilon = rmvl(n,mu = rep(0,q), Sigma = Sigma0)
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}
  

#Multivariate log-Normal
Gen_data9 = function(n = 100, beta, Sigma0){
  q = length(beta)
  #X = rnorm(n,5,1)
  #X = rcauchy(n)
  X = rnorm(n)
  epsilon = exp(rmvn(n,mu = rep(0,q), Sigma = Sigma0))
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}

Gen_data9X = function(n = 100, beta, Sigma0){
  q = length(beta)
  #X = rnorm(n,5,1)
  #X = rcauchy(n)
  X = rnorm(n)
  lep = exp(rnorm(n, 0, Sigma0))
  for (j in 2:q){
    lep = cbind(lep, exp(rnorm(n, 0, Sigma0)))
  }
  #epsilon = exp(rmvn(n,mu = rep(0,q), Sigma = Sigma0))
  Y = X%*%t(beta) + lep
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}
#Gen_data9X(n=100, beta = c(1,0,1), Sigma0 = 2) 


#g-h distribution
library(gk)

Gen_dataTgh = function(n = 100, beta, A, B, g=1, h=0.2){ #as the reviewer requested
  q = length(beta)
  #X = rnorm(n,5,1)
  #X = rcauchy(n)
  X = rnorm(n)
  epsilon = matrix(rep(0, n*q), ncol = q)
  for(j in 1:q){
    epsilon[,j] = rgh(n, A, B, g=1, h=0.2, c = 0.8, type = "tukey")
  }
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}

#Gen_dataTgh(n=100, beta=c(1,0,0), A=0,B=1)

Gen_dataGgh = function(n = 100, beta, A, B, g=1, h=0.2){ #as the reviewer requested
  q = length(beta)
  #X = rnorm(n,5,1)
  #X = rcauchy(n)
  X = rnorm(n)
  epsilon = matrix(rep(0, n*q), ncol = q)
  for(j in 1:q){
    epsilon[,j] = rgh(n, A, B, g=1, h=0.2, c = 0.8, type = "generalised")
  }
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}

