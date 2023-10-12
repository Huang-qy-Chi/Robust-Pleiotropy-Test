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
S0 = diag(10)
S1 = 0.2*matrix(rep(1,9),ncol = 3)+0.8*diag(3)
S2 = 0.5*matrix(rep(1,9),ncol = 3)+0.5*diag(3)
S3 = 0.8*matrix(rep(1,9),ncol = 3)+0.2*diag(3)

cl <- makeCluster(18)
registerDoParallel(cl)
RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
  #result = Gen_data(n = 1000, beta = c(1,0,rep(0,1)), Sigma0 = 2*diag(3))
  #result = Gen_data2(n = 300, beta = c(0.25,0.25,0), S = diag(3),COR = S1)
  #result = Gen_data3(n = 500, beta = c(0.25,0.25,0),mu = rep(10,3), Sigma1 = S1,Sigma2 = 10*diag(3))
  #result = Gen_data5(n = 500, beta = c(1,0,0), Sigma0 = S1)
  result = Gen_data9(n = 500, beta = c(1,0,0), Sigma0 = S1)
  y = result$Y
  x = result$X
  ini = M_LTS2(y,x)
  B0 = ini$beta
  Sigma0 = ini$Sigma
  #MM_estimate1X(y,x,B0,Sigma0)$s
  Me = MM_estimate1X(y,x,B0,Sigma0)
  #Me = MM_estimate1X(y,x,B0,Sigma0)
  BM = Me$Beta
  SigmaM = Me$Sigma
  s=Me$s
  EHT_PFCX(y,x,BM,SigmaM,s,alpha = 0.05)$index
}
mean(RR)
stopImplicitCluster()
stopCluster(cl)
beep(1)




help(package = "LaplacesDemon")













#main:q=3
cl <- makeCluster(18)
registerDoParallel(cl)
RS = matrix(rep(1,12),ncol = 4)
n = c(200,300,500)
corr = c(0.2,0.5,0.8)
alpha = c(0.05,0.01)
Beta = matrix(c(rep(0,3),c(1,rep(0,2)),c(rep(0.25,2),0),rep(0.25,3)),ncol = 4)
sigma = c(1, 2)
#Normal
for(p in 1:2){ #alpha
  for(r in 1:3){ #COrr
    S = corr[r]*matrix(rep(1,9),ncol = 3)+(1-corr[r])*diag(3)
    for(m in 1:3){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          y = result$Y
          x = result$X
          ini = M_LTS2(y,x)
          B0 = ini$beta
          Sigma0 = ini$Sigma
          Me = MM_estimate1X(y,x,B0,Sigma0)
          BM = Me$Beta
          SigmaM = Me$Sigma
          EHT_PFCX(y,x,BM,SigmaM,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500")
    colnames(RS) = c("H00","H01","H10","H11")  
    write.table(RS,sep=",",file =paste("C:/Users/txw/Desktop/EHT/A4/EHTC-3-",100*alpha[p],"/Normal3-",10*corr[r],".csv",sep=""),row.names = T)
    #beep(1)
    print("wait Normal.")
    Sys.sleep(20)
  }
}


#Heteroscedasticity
for(p in 1:2){ #alpha
  for(r in 1:3){ #COrr
    S = corr[r]*matrix(rep(1,9),ncol = 3)+(1-corr[r])*diag(3)
    for(m in 1:3){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          y = result$Y
          x = result$X
          ini = M_LTS2(y,x)
          B0 = ini$beta
          Sigma0 = ini$Sigma
          Me = MM_estimate1X(y,x,B0,Sigma0)
          BM = Me$Beta
          SigmaM = Me$Sigma
          EHT_PFCX(y,x,BM,SigmaM,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500")
    colnames(RS) = c("H00","H01","H10","H11")  
    write.table(RS,sep=",",file =paste("C:/Users/txw/Desktop/EHT/A4/EHTC-3-",100*alpha[p],"/heNormal3-",10*corr[r],".csv",sep=""),row.names = T)
    #beep(1)
    print("wait Hetero.")
    Sys.sleep(20)
  }
}


#Mix Normal
for(p in 1:2){ #alpha
  for(r in 1:3){ #COrr
    S = corr[r]*matrix(rep(1,9),ncol = 3)+(1-corr[r])*diag(3)
    for(m in 1:3){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          result = Gen_data3(n = n[m], beta = Beta[,u],mu = rep(3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          y = result$Y
          x = result$X
          ini = M_LTS2(y,x)
          B0 = ini$beta
          Sigma0 = ini$Sigma
          Me = MM_estimate1X(y,x,B0,Sigma0)
          BM = Me$Beta
          SigmaM = Me$Sigma
          EHT_PFCX(y,x,BM,SigmaM,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500")
    colnames(RS) = c("H00","H01","H10","H11")  
    write.table(RS,sep=",",file =paste("C:/Users/txw/Desktop/EHT/A4/EHTC-3-",100*alpha[p],"/mxNormal3-",10*corr[r],".csv",sep=""),row.names = T)
    #beep(1)
    print("wait Mix.")
    Sys.sleep(20)
  }
}


#Cauchy
for(p in 1:2){ #alpha
  for(r in 1:3){
    S = corr[r]*matrix(rep(1,9),ncol = 3)+(1-corr[r])*diag(3)
    for(m in 1:3){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          result = Gen_data2(n = n[m], beta = Beta[,u], S = diag(3),COR = S)
          y = result$Y
          x = result$X
          ini = M_LTS2(y,x)
          B0 = ini$beta
          Sigma0 = ini$Sigma
          Me = MM_estimate1X(y,x,B0,Sigma0)
          BM = Me$Beta
          SigmaM = Me$Sigma
          EHT_PFCX(y,x,BM,SigmaM,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
  rownames(RS) = c("n=200","n=300","n=500")
  colnames(RS) = c("H00","H01","H10","H11")  
  write.table(RS,sep=",",file =paste("C:/Users/txw/Desktop/EHT/A4/EHTC-3-",100*alpha[p],"/Cauchy3-",10*corr[r],".csv",sep=""),row.names = T)
  #beep(1)
  print("wait Cauchy.")
  Sys.sleep(20)
  }
}


#Log-Normal
for(p in 1:2){ #alpha
  for(r in 1:3){
    S = corr[r]*matrix(rep(1,9),ncol = 3)+(1-corr[r])*diag(3)
    for(m in 1:3){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          result = Gen_data9(n = n[m], beta = Beta[,u], Sigma0 = S)
          y = result$Y
          x = result$X
          ini = M_LTS2(y,x)
          B0 = ini$beta
          Sigma0 = ini$Sigma
          Me = MM_estimate1X(y,x,B0,Sigma0)
          BM = Me$Beta
          SigmaM = Me$Sigma
          EHT_PFCX(y,x,BM,SigmaM,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500")
    colnames(RS) = c("H00","H01","H10","H11")  
    write.table(RS,sep=",",file =paste("C:/Users/txw/Desktop/EHT/A4/EHTC-3-",100*alpha[p],"/Lognormal3-",10*corr[r],".csv",sep=""),row.names = T)
    #beep(1)
    print("wait Lognormal.")
    Sys.sleep(20)
  }
}


for(p in 1:2){ #alpha
  for(r in 1:2){#sigma
    #S = corr[r]*matrix(rep(1,9),ncol = 3)+(1-corr[r])*diag(3)
    for(m in 1:3){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          result = Gen_data9X(n = n[m], beta = Beta[,u], Sigma0 = sigma[r])
          y = result$Y
          x = result$X
          ini = M_LTS2(y,x)
          B0 = ini$beta
          Sigma0 = ini$Sigma
          Me = MM_estimate1X(y,x,B0,Sigma0)
          BM = Me$Beta
          SigmaM = Me$Sigma
          EHT_PFCX(y,x,BM,SigmaM,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500")
    colnames(RS) = c("H00","H01","H10","H11")  
    write.table(RS,sep=",",file =paste("C:/Users/txw/Desktop/EHT/A4/lognormal-3-",100*alpha[p],"/Lognormal3-",sigma[r],".csv",sep=""),row.names = T)
    #beep(1)
    print("wait LognormalX.")
    Sys.sleep(20)
  }
}


#g-and-h distribution
library(gk)
for(p in 1:2){ #alpha
  for(r in 1:2){#sigma
    #S = corr[r]*matrix(rep(1,9),ncol = 3)+(1-corr[r])*diag(3)
    for(m in 1:3){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = c("LaplacesDemon","gk"))%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          result = Gen_dataTgh(n = n[m], beta = Beta[,u], A = 0, B = sigma[r])
          y = result$Y
          x = result$X
          ini = M_LTS2(y,x)
          B0 = ini$beta
          Sigma0 = ini$Sigma
          Me = MM_estimate1X(y,x,B0,Sigma0)
          BM = Me$Beta
          SigmaM = Me$Sigma
          EHT_PFCX(y,x,BM,SigmaM,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500")
    colnames(RS) = c("H00","H01","H10","H11")  
    write.table(RS,sep=",",file =paste("C:/Users/txw/Desktop/EHT/A4/lognormal-3-",100*alpha[p],"/Tgh3-",sigma[r],".csv",sep=""),row.names = T)
    #beep(1)
    print("wait Tgh.")
    Sys.sleep(20)
  }
}


for(p in 1:2){ #alpha
  for(r in 1:2){#sigma
    #S = corr[r]*matrix(rep(1,9),ncol = 3)+(1-corr[r])*diag(3)
    for(m in 1:3){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = c("LaplacesDemon","gk"))%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          #result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          result = Gen_dataGgh(n = n[m], beta = Beta[,u], A = 0, B = sigma[r])
          y = result$Y
          x = result$X
          ini = M_LTS2(y,x)
          B0 = ini$beta
          Sigma0 = ini$Sigma
          Me = MM_estimate1X(y,x,B0,Sigma0)
          BM = Me$Beta
          SigmaM = Me$Sigma
          EHT_PFCX(y,x,BM,SigmaM,alpha = alpha[p])$index
        }
        RS[m,u]=mean(RR)
      }
    }
    rownames(RS) = c("n=200","n=300","n=500")
    colnames(RS) = c("H00","H01","H10","H11")  
    write.table(RS,sep=",",file =paste("C:/Users/txw/Desktop/EHT/A4/lognormal-3-",100*alpha[p],"/Ggh3-",sigma[r],".csv",sep=""),row.names = T)
    #beep(1)
    print("wait Ggh.")
    Sys.sleep(20)
  }
}






stopImplicitCluster()
stopCluster(cl)
print("Mission accomplish!")
#beep(3)
#system("shutdown /h/ t 6/f")
#shell.exec("C:\\WINDOWS\\SYSTEM32\\shutdown.exe")
#shell("shutdown /h/ t 6/f",translate = T)


























cl <- makeCluster(7)
registerDoParallel(cl)
RS = matrix(rep(0,16),ncol = 4)
n = c(200,300,500,1000)
Beta = matrix(c(rep(0,10),c(1,rep(0,9)),c(rep(1,5),rep(0,5)),rep(1,10)),ncol = 4)
for(m in 1:4){
  for(u in 1:4){
    RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
      #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S3,Sigma2 = 3*diag(3))
      #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S3)
      result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S1)
      #result = Gen_data2(n = n[m], beta = Beta[,u], S = S0, COR =S2)
      y = result$Y
      x = result$X
      ini = M_LTS2(y,x)
      B0 = ini$beta
      Sigma0 = ini$Sigma
      Me = MM_estimate1(y,x,B0,Sigma0)
      BM = Me$Beta
      SigmaM = Me$Sigma
      EHT_PFC(y,x,BM,SigmaM,alpha = 0.05)$index
    }
    RS[m,u]=mean(RR)
  }
}
rownames(RS) = c("n=200","n=300","n=500","n=1000")
colnames(RS) = c("H00","Hk0","H10","H11")  
stopImplicitCluster()
stopCluster(cl)
write.table(RS,sep=",",file ="C:/Users/Administrator/Desktop/EHT/EHTC-10-5/Normal2.csv",row.names = T)
beep(1)







