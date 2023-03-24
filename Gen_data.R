# Function of Data Genetate
library(MASS)
library(LaplacesDemon)
# Case1: Vec(e^T)~N(0,In\otimes Sigma)
Gen_data1 = function(n = 100, beta, Sigma0){
  q = ncol(Sigma0)
  In = diag(rep(1,n))
  vecVar = kronecker(In,Sigma0)
  VE = MASS::mvrnorm(n*q, mu = rep(0,n*q), Sigma = vecVar)#按行拉直后重组矩???
  E = matrix(VE,nrow = n, ncol = q, byrow = T) #矩阵正态扰???
  X = matrix(rnorm(n),ncol = 1)
  Y = X[1,]%*%beta + E[1,]
  for(i in 2:n){
    Y = c(Y, X[i,]%*%beta + E[i,]) #按行拉直的Y
  }
  Y = matrix(Y, nrow = n, ncol = q, byrow = T)
  Y_string = sapply(1:q, function(u){return(paste("Y", u, sep = ""))})
  colnames(Y) <- Y_string
  colnames(X) <- c("X")
  return(list(Y = Y, X = X))
}


library(LaplacesDemon)
Gen_data = function(n = 100, beta, Sigma0){
  q = length(beta)
  X = rnorm(n,0,1)
  epsilon = rmvn(n, mu = rep(0,q), Sigma = Sigma0)
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}


Gen_data2 = function(n = 100, beta, S, COR){
  library(LaplacesDemon)
  q = length(beta)
  X = rnorm(n)
  epsiolon = rmvc(n, mu = rep(0,q), S = S)%*%Msq(COR,0.5)
  Y = X%*%t(beta) + epsiolon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}

Gen_data3 = function(n = 100, beta, mu, Sigma1, Sigma2,ratio = 0.2){
  library(LaplacesDemon)
  q = as.numeric(ncol(Sigma1))
  X = rnorm(n)
  w =runif(1)
  if(w<1-ratio){
    epsilon = rmvn(n, mu = rep(0,q), Sigma = Sigma1)
  }else{
    epsilon = rmvn(n, mu = mu, Sigma = Sigma2)
  }
    Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}

#Question
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
Gen_data5 = function(n = 100, beta, Sigma0){
  q = as.numeric(ncol(Sigma0))
  X = rnorm(n)
  epsilon = (diag(n)+1*diag(X))%*%rmvn(n, mu = rep(0,q), Sigma = Sigma0)
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}




Gen_data6 = function(n = 100, beta, Sigma0){
  q = length(beta)
  #X = rnorm(n,5,1)
  #X = rcauchy(n)
  X = rexp(n)
  epsilon = rmvn(n, mu = rep(0,q), Sigma = Sigma0)
  Y = X%*%t(beta) + epsilon
  return(list(Y = matrix(Y,ncol = q), X = matrix(X,ncol = 1)))
}






















