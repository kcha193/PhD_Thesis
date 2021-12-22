


update.Infor.mat =
function(int, trt.mat, blk.proj){

  library(MASS)

  Ab = t(trt.mat) %*% blk.proj %*% trt.mat
  ginv.Ab = ginv(Ab)
  Rep = nrow(trt.mat)/ncol(trt.mat)

  a = int[1]
  b = int[2]

  T1 = trt.mat

  x1 = T1[b,] = trt.mat[a,]
  x2 = T1[a,] = trt.mat[b,]

  if(isTRUE(all.equal(x1,x2))) return(NULL)

  c1 = blk.proj[a,a]
  c2 = blk.proj[b,b]

  d1 = blk.proj[a, -c(a,b)]
  d2 = blk.proj[b, -c(a,b)]

  T2 = trt.mat[-c(a,b),]


  B = (x1 - x2) %*% t(d1 - d2) %*% T2

  #return( list( n.trt.mat = T1,
  #A = (t(trt.mat) %*% blk.proj %*% trt.mat) -
  Aa = Ab - (((c1-c2) * (x1%*% t(x1) -x2%*%t(x2))) + B + t(B))


  e.va = eigen(Aa)$va
  Aa.ave.eff = (length(e.va[-which(e.va<1e-7)])/sum(1/e.va[-which(e.va<1e-7)]))/Rep

  e.va = eigen(Ab)$va
  Ab.ave.eff = (length(e.va[-which(e.va<1e-7)])/sum(1/e.va[-which(e.va<1e-7)]))/Rep


  if(Ab.ave.eff>=  Aa.ave.eff){
    return (NULL)
  }else{
    return(list(AbAveEff = Ab.ave.eff,
                newTrtMat = T1,
                AaAveEff = Aa.ave.eff))

  }
}
