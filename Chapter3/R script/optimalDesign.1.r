sourceDir <- function(path, trace = TRUE, ...) {
  library(MASS)
  library(lattice)
  library(Rcpp)
  library(inline)
  library(RcppArmadillo)
  library(compiler)

  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
  }
}



test = 
function(X.trt, blk.proj, C.trt.mat = diag(ncol(X.trt)), Rep){
  if( length(dim(X.trt))==3) X.trt = X.trt[,,1]
  
  info.mat = C.trt.mat %*% t(X.trt) %*% blk.proj %*% X.trt %*% C.trt.mat
  trace = tr(info.mat)
  e.va = eigen(info.mat)$va
  e.vec =  eigen(info.mat)$vec
  can.eff = e.va[-which(e.va<1e-7)]/Rep
  list( trace = trace,
        nCan = length(can.eff),
        can.eff = can.eff,
        ave.eff =  1/mean(1/can.eff),
        e.vec = e.vec )
}

sourceDir(path = "C:/Users/kcha193/My Dropbox/R functions for two-phase experiments")
sourceDir(path = "C:/Users/Kevin/Documents/My Dropbox/R functions for two-phase experiments")

#Setup the design

#Parameter's of block
n = 100
nBlk = 20
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
nTrt = 20
Rep = (nBlk*nPlot)/nTrt

trt.des = as.factor(rep(1:nTrt, Rep))
X.trt = matrix(0, ncol = nTrt, nrow= nTrt*Rep)
X.trt[cbind(1:( nTrt*Rep), trt.des)] = 1

test(X.trt, blk.proj, Rep) 

design.df$temp = as.factor(apply(X.trt, 1, function(x) which(as.logical(x))))

summary.aov.onePhase(design.df,  blk.str = "Block" , trt.str = "temp")

#Randomised method:
Rprof("example.out")
ccp.r.opt.X.trt = cpp.r.optimised(X.trt, blk.proj, Rep, 10000000)
Rprof(NULL)
summaryRprof("example.out") 
test(ccp.r.opt.X.trt$X.trt, blk.proj, Rep) 

design.df$temp = as.factor(apply(ccp.r.opt.X.trt$X.trt, 1, function(x) 
which(as.logical(x))))

summary.aov.onePhase(design.df,  blk.str = "Block" , trt.str = "temp")


cmpfun.r.optimised = cmpfun(r.optimised)
Rprof("example.out")
r.opt.X.trt = cmpfun.r.optimised(X.trt, blk.proj)
Rprof(NULL)
summaryRprof("example.out") 

test(r.opt.X.trt, blk.proj, Rep = Rep) 

design.df$temp = as.factor(apply(r.opt.X.trt, 1, function(x) which(as.logical(x))))

summary.aov.onePhase(design.df,  blk.str = "Block" , trt.str = "temp")

#Swap method:
comfun.s.optimised = cmpfun(s.optimised)
#t1 = proc.time()
Rprof("example.out")
s.opt.X.trt = s.optimised(X.trt, blk.proj, nIter = 10)
#proc.time()-t1
Rprof(NULL)
summaryRprof("example.out") 

X.trt = s.opt.X.trt[[1]]

test(X.trt, blk.proj, Rep) 

design.df$temp = as.factor(apply(X.trt, 1, function(x) which(as.logical(x))))

summary.aov.onePhase(design.df,  blk.str = "Block" , trt.str = "temp")


#Swap method:
#t1 = proc.time()
Rprof("example.out")
cpp.s.opt.X.trt = cpp.s.optimised(X.trt, blk.proj, nIter = 10)
#proc.time()-t1
Rprof(NULL)
summaryRprof("example.out") 

X.trt = cpp.s.opt.X.trt[,,1]

test(X.trt, blk.proj, Rep) 

design.df$temp = as.factor(apply(X.trt, 1, function(x) which(as.logical(x))))

summary.aov.onePhase(design.df,  blk.str = "Block" , trt.str = "temp")

#Swap method:
#t1 = proc.time()
Rprof("example.out")
cpp.new.s.opt.X.trt = cpp.new.s.optimised(X.trt, blk.proj)
#proc.time()-t1
Rprof(NULL)
summaryRprof("example.out") 

test(cpp.new.s.opt.X.trt$X.trt, blk.proj, Rep = Rep) 

design.df$temp = as.factor(apply(cpp.new.s.opt.X.trt$best.X.trt, 1, function(x) which(as.logical(x))))

summary.aov.onePhase(design.df,  blk.str = "Block" , trt.str = "temp")

################################################################################
#Test on 4-by-3 factorial experiments 
#Parameter's of block
n = 10
nBlk = 2
nPlot = n/nBlk

Zb = matrix(0, ncol = nBlk, nrow = n)
Z.des = rep(1:nBlk,  each = nPlot)
Zb[cbind(1:n, Z.des)] <- 1

design.df = data.frame(Block = as.factor(Z.des))

Pb = projMat(Zb)

betBlock = Pb - mK(n)
withBlock = mI(n) - Pb                  

blk.proj = withBlock


#Parameter's of treatment
nTrt1 = 2
nTrt2 = 5    
#nTrt3 = 3    

nTrt = nTrt1 * nTrt2 #* nTrt3

Rep = n/nTrt

trt.des = as.factor(rep(1:nTrt, Rep))
X.trt = matrix(0, ncol = nTrt, nrow= nTrt*Rep)
X.trt[cbind(1:( nTrt*Rep), trt.des)] = 1

C1 = (mI(nTrt1) - mK(nTrt1))  %x%  mK(nTrt2)              %x%  mK(nTrt3)
C2 =  mK(nTrt1)               %x% (mI(nTrt2) - mK(nTrt2)) %x%  mK(nTrt3)
C3 =  mK(nTrt1)               %x%  mK(nTrt2)              %x%  (mI(nTrt3) - mK(nTrt3))

C4 =  mK(nTrt1)               %x% (mI(nTrt2) - mK(nTrt2)) %x% (mI(nTrt3) - mK(nTrt3))
C5 = (mI(nTrt1) - mK(nTrt1))  %x% (mI(nTrt2) - mK(nTrt2)) %x%  mK(nTrt3)
C6 = (mI(nTrt1) - mK(nTrt1))  %x%  mK(nTrt2)              %x% (mI(nTrt3) - mK(nTrt3))

C7 = (mI(nTrt1) - mK(nTrt1))  %x% (mI(nTrt2) - mK(nTrt2)) %x%  (mI(nTrt3) - mK(nTrt3))

C.trt.mat = rbind(C5, C1, C2)


design.df$Trt1 = factor(LETTERS[rep(1:nTrt1, n/nTrt1)])
design.df$Trt2 = factor(c("a", "b", "b", "a"))
design.df$Trt3 = factor(rep(1:nTrt3, n/nTrt3))

design.df$"Trt1:Trt2:Trt3" = factor(rep(1:nTrt, Rep))

getReplicationList(design.df, c("Trt1", "Trt2","Trt3", "Trt1:Trt2", "Trt1:Trt3", "Trt2:Trt3"))

#Swap method:
Rprof("example.out")
cpp.new.s.opt.X.trt = cpp.new.fact.s.optimised(X.trt, blk.proj, C.trt.mat)
Rprof(NULL)
summaryRprof("example.out") 

test(cpp.new.s.opt.X.trt$best.X.trt, blk.proj, C7, Rep)
test(cpp.new.s.opt.X.trt$best.X.trt, blk.proj, C6, Rep)
test(cpp.new.s.opt.X.trt$best.X.trt, blk.proj, C5, Rep)
test(cpp.new.s.opt.X.trt$best.X.trt, blk.proj, C4, Rep)
test(cpp.new.s.opt.X.trt$best.X.trt, blk.proj, C3, Rep)
test(cpp.new.s.opt.X.trt$best.X.trt, blk.proj, C2, Rep)
test(cpp.new.s.opt.X.trt$best.X.trt, blk.proj, C1, Rep)

colnames(cpp.new.s.opt.X.trt$best.X.trt) = sort(levels(interaction(LETTERS[1:nTrt1],letters[1:nTrt2],1:nTrt3)))

trt.des = apply(cpp.new.s.opt.X.trt$best.X.trt, 1,  function(x) 
colnames(cpp.new.s.opt.X.trt$best.X.trt)[which(as.logical(x))])
 trt.des = as.data.frame(t(sapply(strsplit(trt.des, "\\."), function(x) x)))

colnames(trt.des) = c("Trt1", "Trt2", "Trt3")
 
design.df = data.frame(Block = as.factor(Z.des))

design.df = cbind(design.df, trt.des)

summary.aov.onePhase(design.df,  blk.str = "Block" , trt.str = "Trt1 + Trt2 + Trt3")


cpp.new.s.opt.X.trt = cpp.new.s.optimised(X.trt, blk.proj, C2)

test(cpp.new.s.opt.X.trt$X.trt, blk.proj, C3, Rep)
test(cpp.new.s.opt.X.trt$X.trt, blk.proj, C2, Rep)
test(cpp.new.s.opt.X.trt$X.trt, blk.proj, C1, Rep)

################################################################################
################################################################################
#CRD on 

nAni = 6
nTrt = 2
phase1DesignEX1 <- local({
  Ani = as.factor(LETTERS[1:nAni])
  Trt = as.factor(letters[1:nTrt])
  data.frame(Ani,Trt)
})
phase1DesignEX1

summary.aov.onePhase(phase1DesignEX1,  blk.str = "Ani",  trt.str = "Trt")

#Parameter's of block
n = 2
nPlot = 4
nBlk = n/nPlot

Zb = matrix(0, ncol = nBlk, nrow = n)
Z.des = rep(1:nBlk,  each = nPlot)
Zb[cbind(1:n, Z.des)] <- 1

Zt = matrix(0, ncol = nPlot, nrow = n)
tag.des = rep(1:nPlot,  time = nBlk)
Zt[cbind(1:n, tag.des)] <- 1

Pb = projMat(Zb)
Pb1 = projMat(Zt)

betRun = Pb - mK(n)
betTag = Pb1 - mK(n)
withBlock = (mI(n) - Pb) # %*% (mI(n) - Pb1)

blk.proj = withBlock


#Parameter's of treatment
nTrt = nrow(phase1DesignEX1)

Rep = n/nTrt

trt.des = as.factor(rep(1:nTrt, Rep))
X.trt = matrix(0, ncol = nTrt, nrow= nTrt*Rep)
X.trt[cbind(1:( nTrt*Rep), trt.des)] = 1

xxxx = cpp.new.fact.s.optimised(X.trt, blk.proj)

test(X.trt = xxxx$best.X.trt, blk.proj, Rep = Rep) 

colnames(xxxx$best.X.trt) = sort(levels(interaction(LETTERS[1:nTrt])))

trt.des = apply(xxxx$best.X.trt, 1,  function(x) 
colnames(xxxx$best.X.trt)[which(as.logical(x))])
 trt.des = as.data.frame(t(sapply(strsplit(trt.des, "\\."), function(x) x)))

design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

design.df = cbind(design.df, 
            Ani = t(trt.des), 
            Trt = phase1DesignEX1[match(as.character(t(trt.des)), 
                          as.character(phase1DesignEX1$Ani)),]$Trt)

design.df

summary.aov.onePhase(design.df,  blk.str = "Run + Tag" , trt.str = "Ani")

summary.aov.twoPhase(design.df,  blk.str2 = "Run" , blk.str1 = "Ani", trt.str = "Trt + Tag")

summary.aov.twoPhase(design.df,  blk.str2 = "Run + Tag" , blk.str1 = "Ani", trt.str = "Trt")

################################################################################
################################################################################
#RCD on 

#Phase 2 Block structure 
n = 12
nPlot = 4
nBlk = n/nPlot

Zb = matrix(0, ncol = nBlk, nrow = n)
Z.des = rep(1:nBlk,  each = nPlot)
Zb[cbind(1:n, Z.des)] <- 1

design.df = data.frame(Block = as.factor(Z.des))

Zt = matrix(0, ncol = nPlot, nrow = n)
tag.des = rep(1:nPlot,  time = nBlk)
Zt[cbind(1:n, tag.des)] <- 1

Pb = projMat(Zb)
Pb1 = projMat(Zt)

betRun = Pb - mK(n)
betTag = Pb1 - mK(n)     
withBlock = (mI(n) - Pb) %*% (mI(n) - Pb1)                  

blk.proj = withBlock


#Parameter's of Phase 1 BLock structure 
nCg = 2
nAni = 3    

nB1 = nCg * nAni

Rep = n/nB1

b1.des = as.factor(rep(1:nB1, Rep))
B1.trt = matrix(0, ncol = nB1, nrow= nB1*Rep)
B1.trt[cbind(1:( nB1*Rep), b1.des)] = 1

C1 = (mI(nCg) - mK(nCg))  %x%  mK(nAni)             
C2 =  mI(nCg)             %x% (mI(nAni) - mK(nAni))  

C.trt.mat = rbind(C2, C1)

#Swap method:
Rprof("example.out")
cpp.new.s.opt.X.trt = cpp.new.fact.s.optimised(B1.trt, blk.proj, C.trt.mat)
Rprof(NULL)
summaryRprof("example.out") 
cpp.new.s.opt.X.trt
test(cpp.new.s.opt.X.trt$best.X.trt, blk.proj, C2, Rep)

colnames(cpp.new.s.opt.X.trt$best.X.trt) = sort(levels(interaction(LETTERS[1:nCg],letters[1:nAni])))

trt.des = apply(cpp.new.s.opt.X.trt$best.X.trt, 1,  function(x) 
colnames(cpp.new.s.opt.X.trt$best.X.trt)[which(as.logical(x))])
 trt.des = as.data.frame(t(sapply(strsplit(trt.des, "\\."), function(x) x)))

 colnames(trt.des) = c("Cage", "Ani")
 
design.df = data.frame(Block = as.factor(Z.des))

design.df = cbind(design.df, trt.des)

summary.aov.onePhase(design.df,  blk.str = "Block" , trt.str = "Cage/Ani")

summary.aov.onePhase(design.df,  blk.str = "Block" , trt.str = "Cage:Ani")

with(design.df, table(Block, Cage, Ani))

test(cpp.new.s.opt.X.trt$best.X.trt, blk.proj, C2, Rep)
test(cpp.new.s.opt.X.trt$best.X.trt, blk.proj, C1, Rep)

################################################################################
# A block design for 7 treatments in 7 blocks of size 3. Note how withinData 
# is recycled to fill out the blocksize requirements.

nBlk = 6
nPlot = 4
n = nBlk * nPlot 

Zb = matrix(0, ncol = nBlk, nrow = n)
Z.des = rep(1:nBlk,  nPlot)
Zb[cbind(1:n, Z.des)]<-1

design.df = data.frame(Block = as.factor(Z.des))

Pb = projMat(Zb)

betBlock = Pb - mK(n)
withBlock = mI(n) - Pb                  

blk.proj = withBlock

#Parameter's of treatment
nTrt = 6
Rep = (nBlk*nPlot)/nTrt

trt.des = as.factor(rep(1:nTrt, Rep))
X.trt = matrix(0, ncol = nTrt, nrow= nTrt*Rep)
X.trt[cbind(1:( nTrt*Rep), trt.des)] = 1


#Cutting algorithm
S = t(Zb) %*% X.trt
sj = matrix(colSums(S), nrow = nBlk, ncol = nTrt, byrow =TRUE)
newS = S - (nBlk/n)*sj

sum(newS^2)


ginv(t(X.trt) %*% withBlock %*% X.trt)

xc  = X.trt - 1/ncol(X.trt)

Xc = cbind(Zb, xc)

solve(t(Xc) %*% Xc)

det(t(Xc) %*% Xc)

det(t(Zb) %*% Zb)

det(ginv(t(X.trt) %*% withBlock %*% X.trt))

library(AlgDesign)

BIB<-optBlock(~.,withinData=factor(1:nTrt),blocksizes=rep(nPlot,nBlk), criterion="D")

design.df$temp =  BIB$design$X1
summary.aov.onePhase(design.df,  blk.str = "Run" , trt.str = "temp")

BIB<-optBlock(~.,withinData=factor(1:15), blocksizes=rep(5,42), criterion="OB")
design.df$temp =  BIB$design$X1
summary.aov.onePhase(design.df,  blk.str = "Block" , trt.str = "temp")


crossprod(table(c(rep(1:42, rep(5,42))),BIB$design[,1]))

nBlk = 8
nPlot = 8
n = nBlk * nPlot


Zb = matrix(0, ncol = nBlk, nrow = n)
Z.des = rep(1:nBlk,  each = nPlot)
Zb[cbind(1:n, Z.des)]<-1

design.df = data.frame(Block = as.factor(Z.des))


dat = data.frame(Cag = factor(rep(1:8, each = 8)),
                Ani = factor(rep(1:2, each = 4)),
                Sam = factor(rep(1:2, each = 2)),
                Sub = factor(1:2))


des = optBlock(~.,withinData=dat,blocksizes=rep(8,8), criterion="OB")

design.df=  cbind(design.df, des$design)
summary.aov.onePhase(design.df,  blk.str = "Block" , trt.str = "Cag/Ani/Sam/Sub")


Tag = gen.factorial(8,1,factor =1, varName = "Tag")

final<-optBlock(~.,withinData=Tag,wholeBlockData=des$Blocks, blocksizes=rep(8,8), criterion="D")

des<-optBlock(~.,withinData=Tag, blocksizes=rep(8,8), criterion="D")

final<-optBlock(~.,withinData=dat, wholeBlockData=des$Blocks, blocksizes=rep(1,8), criterion="D")
