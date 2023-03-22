# Function to derive a matrix's n square:based on eigen-decomposition
Msq = function(x, a){
  m = ncol(x)
  ei = eigen(x)$values
  eigen0 = diag(ei)
  Vector = eigen(x)$vectors
  if(a<1){
    for(i in 1:m){
      if(ei[i]<0){
        print("Not positive definite!")
        break
      }
    }
  }
  ei1 = ei^(a)
  eigen1 = diag(ei1)
  ed = (Vector)%*%eigen1%*%t(Vector)
  return(ed)
}
