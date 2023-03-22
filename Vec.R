#把以矩阵存储的数据拉直,应该注意是行数据还是列式列数据,函数内部默认为行数据
Vec = function(x, byrow=T, is.na=F){#其中x为矩阵,byrow为矩阵数据是按行还是按列拉直                                       #is.na表示是否要处理原数据中的确实值,其值赋为T时，表示删除缺失值
  result <- vector()
  if(byrow==T){   #按行拉直
    for(i in 1:nrow(x)){
      result<-c(result,as.numeric(x[i,]))
    }
  }
  if(byrow==F){   #按列拉直
    for(i in 1:ncol(x)){
      result<-c(result,as.numeric(x[,i]))
    }
  }
  if(is.na==T){  
    return(na.omit(result))
  }else{
    return(result)
    }
}
