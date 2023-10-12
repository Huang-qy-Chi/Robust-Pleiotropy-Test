#Main function
rm(list = ls())
library(LaplacesDemon)
library(foreach)
library(doParallel)
library(beepr)
setwd("C:/Users/txw/Desktop/EHT")
source("Gen_data.R")
source("Vec.R")
source("M_LTSfast.R")
source("MM_estimatefast.R")
source("Msq.R")
source("EHT_Pfast.R")



Sys.sleep(120)
source("EHT_MLH.R")

#main:q=10
cl <- makeCluster(6)
registerDoParallel(cl)
RS = matrix(rep(0,16),ncol = 4)
n = c(200,300,500,1000)
corr = c(0.2,0.5,0.8)
alpha = c(0.05,0.01)
Beta = matrix(c(rep(0,10),c(1,rep(0,9)),c(rep(0.25,2),rep(0,8)),rep(0.25,10)),ncol = 4)

#Normal
for(p in 1:2){ #alpha
  for(r in 1:3){ #COrr
    S = corr[r]*matrix(rep(1,100),ncol = 10)+(1-corr[r])*diag(10)
    for(m in 1:4){ #sample size
      for(u in 1:2){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          y = result$Y
          x = result$X
          EHT_MLH(y,x,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500","n=1000")
    colnames(RS) = c("H00","H01","H10","H11")  
    write.table(RS,sep=",",file =paste("C:/Users/Administrator/Desktop/EHT/A1/EHTM-10-",100*alpha[p],"/Normal",10*corr[r],".csv",sep=""),row.names = T)
    beep(1)
    Sys.sleep(20)
  }
}


#Heteroscedasticity
for(p in 1:2){ #alpha
  for(r in 1:3){ #COrr
    S = corr[r]*matrix(rep(1,100),ncol = 10)+(1-corr[r])*diag(10)
    for(m in 1:4){ #sample size
      for(u in 1:2){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          y = result$Y
          x = result$X
          EHT_MLH(y,x,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500","n=1000")
    colnames(RS) = c("H10","H11")  
    write.table(RS,sep=",",file =paste("C:/Users/Administrator/Desktop/EHT/A1/EHTM-10-",100*alpha[p],"/heNormal",10*corr[r],".csv",sep=""),row.names = T)
    beep(1)
    Sys.sleep(20)
  }
}


#Mix Normal
for(p in 1:2){ #alpha
  for(r in 1:3){ #COrr
    S = corr[r]*matrix(rep(1,100),ncol = 10)+(1-corr[r])*diag(10)
    for(m in 1:4){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          result = Gen_data3(n = n[m], beta = Beta[,u],mu = rep(10,10), Sigma1 = S,Sigma2 = 10*diag(10))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          y = result$Y
          x = result$X
          EHT_MLH(y,x,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500","n=1000")
    colnames(RS) = c("H00","H01","H10","H11")   
    write.table(RS,sep=",",file =paste("C:/Users/Administrator/Desktop/EHT/A1/EHTM-10-",100*alpha[p],"/mxNormal",10*corr[r],".csv",sep=""),row.names = T)
    beep(1)
    Sys.sleep(20)
  }
}

#Cauchy
for(p in 1:2){ #alpha
  for(r in 1:3){
    S = corr[r]*matrix(rep(1,100),ncol = 10)+(1-corr[r])*diag(10)
    for(m in 1:4){ #sample size
      for(u in 1:2){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          result = Gen_data2(n = n[m], beta = Beta[,u], S = diag(10),COR = S)
          y = result$Y
          x = result$X
          EHT_MLH(y,x,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500","n=1000")
    colnames(RS) = c("H10","H11")  
    write.table(RS,sep=",",file =paste("C:/Users/Administrator/Desktop/EHT/A1/EHTM-10-",100*alpha[p],"/Cauchy",10*corr[r],".csv",sep=""),row.names = T)
    beep(1)
    Sys.sleep(20)
  }
  
}

stopImplicitCluster()
stopCluster(cl)
beep(3)

Sys.sleep(120)

#main:q=3
cl <- makeCluster(6)
registerDoParallel(cl)
RS = matrix(rep(0,16),ncol = 4)
n = c(200,300,500,1000)
corr = c(0.2,0.5,0.8)
alpha = c(0.05,0.01)
Beta = matrix(c(rep(0,3),c(1,0,0),c(rep(0.25,2),0),rep(0.25,3)),ncol = 4)

#Normal
for(p in 1:2){ #alpha
  for(r in 1:3){ #COrr
    S = corr[r]*matrix(rep(1,9),ncol = 3)+(1-corr[r])*diag(3)
    for(m in 1:4){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          y = result$Y
          x = result$X
          EHT_MLH(y,x,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500","n=1000")
    colnames(RS) = c("H00","H01","H10","H11")
    write.table(RS,sep=",",file =paste("C:/Users/Administrator/Desktop/EHT/A1/EHTM-3-",100*alpha[p],"/Normal",10*corr[r],".csv",sep=""),row.names = T)
    beep(1)
    Sys.sleep(20)
  }
}


#Heteroscedasticity
for(p in 1:2){ #alpha
  for(r in 1:3){ #COrr
    S = corr[r]*matrix(rep(1,9),ncol = 3)+(1-corr[r])*diag(3)
    for(m in 1:4){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          y = result$Y
          x = result$X
          EHT_MLH(y,x,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500","n=1000")
    colnames(RS) = c("H00","H01","H10","H11") 
    write.table(RS,sep=",",file =paste("C:/Users/Administrator/Desktop/EHT/A1/EHTM-3-",100*alpha[p],"/heNormal",10*corr[r],".csv",sep=""),row.names = T)
    beep(1)
    Sys.sleep(20)
  }
}


#Mix Normal
for(p in 1:2){ #alpha
  for(r in 1:3){ #COrr
    S = corr[r]*matrix(rep(1,9),ncol = 3)+(1-corr[r])*diag(3)
    for(m in 1:4){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          result = Gen_data3(n = n[m], beta = Beta[,u],mu = rep(10,3), Sigma1 = S,Sigma2 = 10*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          y = result$Y
          x = result$X
          EHT_MLH(y,x,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500","n=1000")
    colnames(RS) = c("H00","H01","H10","H11") 
    write.table(RS,sep=",",file =paste("C:/Users/Administrator/Desktop/EHT/A1/EHTM-3-",100*alpha[p],"/mxNormal",10*corr[r],".csv",sep=""),row.names = T)
    beep(1)
    Sys.sleep(20)
  }
}

#Cauchy
for(p in 1:2){ #alpha
  for(r in 1:3){
    S = corr[r]*matrix(rep(1,9),ncol = 3)+(1-corr[r])*diag(3)
    for(m in 1:4){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          result = Gen_data2(n = n[m], beta = Beta[,u], S = diag(3),COR = S)
          y = result$Y
          x = result$X
          EHT_MLH(y,x,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500","n=1000")
    colnames(RS) = c("H00","H01","H10","H11")  
    write.table(RS,sep=",",file =paste("C:/Users/Administrator/Desktop/EHT/A1/EHTM-3-",100*alpha[p],"/Cauchy",10*corr[r],".csv",sep=""),row.names = T)
    beep(1)
    Sys.sleep(20)
  }
  
}

stopImplicitCluster()
stopCluster(cl)

Sys.sleep(120)
beep(3)

cl <- makeCluster(6)
registerDoParallel(cl)
RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
  #result = Gen_data(n = 300, beta = c(0.25,0.25,rep(0,8)), Sigma0 = S1)
  #result = Gen_data2(n = 300, beta = c(0.25,0.25,0), S = diag(3),COR = S1)
  #result = Gen_data3(n = 500, beta = c(0.25,0.25,0),mu = rep(10,3), Sigma1 = S1,Sigma2 = 10*diag(3))
  result = Gen_data5(n = 200, beta = c(1,0,0), Sigma0 = S1)
  y = result$Y
  x = result$X
  #ini = M_LTS2(y,x)
  #B0 = ini$beta
  #Sigma0 = ini$Sigma
  #Me = MM_estimate1(y,x,B0,Sigma0)
  #BM = Me$Beta
  #SigmaM = Me$Sigma
  EHT_MLH(y,x,alpha = 0.05)
}
mean(RR)
stopImplicitCluster()
stopCluster(cl)
beep(1)
















