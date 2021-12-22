sourceDir <- function(path, trace = TRUE, ...) {
  library(MASS)
  library(lattice)
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
  }
}



test =
function(X.trt, blk.proj, Rep){
  trace = tr(t(X.trt) %*% blk.proj %*%X.trt)
  e.va = eigen(t(X.trt) %*% blk.proj %*%X.trt)$va

  list(trace = trace,
      con.eff = e.va[-which(e.va<1e-7)]/Rep,
      ave.eff = (length(e.va[-which(e.va<1e-7)])/sum(1/e.va[-which(e.va<1e-7)]))/Rep)
}


sourceDir(path = "C:/Users/kcha193/My Dropbox/R functions for two-phase experiments")


#Parameter's of block
n = 210
nBlk = 42
nPlot = n/nBlk

Zb = matrix(0, ncol = nBlk, nrow = n)
Z.des = rep(1:nBlk,  each = nPlot)
Zb[cbind(1:n, Z.des)]<-1

design.df = data.frame(Block = as.factor(Z.des))

Pb = projMat(Zb)

betBlock = Pb - mK(n)
withBlock = mI(n) - Pb                  

blk.proj = withBlock

#Parameter's of treatment
nTrt = 15
Rep = (nBlk*nPlot)/nTrt

trt.des = as.factor(rep(1:nTrt, Rep))
X.trt = matrix(0, ncol = nTrt, nrow= nTrt*Rep)
X.trt[cbind(1:( nTrt*Rep), trt.des)] = 1


#Construct method
#t1 = proc.time()
c.opt.X.trt = c.optimised(X.trt, blk.proj)
#proc.time()-t1

