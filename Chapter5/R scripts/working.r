# 4/04/2013 11:23:13 a.m.
sourceDir <- function(path, trace = TRUE, ...) {
  library(MASS)
  library(inline)
  library(compiler)
  library(formatR)
  library(Rcpp)
  library(RcppArmadillo)
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
  }
}

sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R")
sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/optimTE/R")

old = proc.time()
design.df = optCRD(nTrt  = 2, bRep  = 3, tRep  = 2, nPlot = 4, iter  = 1000)
proc.time() - old


design.summary.CRD(design.df, FALSE)

sim1 = simDataChisq(design =  design.df, 
            blk.str2 = "Run", blk.str1 = "Ani", 
            trt.str = "Tag + Trt", 
            gamma.R = c(0, 0.25, 1, 4, 100),                    
            gamma.A =c(10^((-8:8)/2)), 
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 2, neg.VC = TRUE, average = FALSE)                        

edfPlot(sim1)

sim2 = simDataChisq(design =  design.df, 
            blk.str2 = "Run", blk.str1 = "Ani", 
            trt.str = "Tag + Trt", 
            gamma.R = c(0, 0.25, 1, 4, 100),                    
            gamma.A =c(10^((-8:8)/2)), 
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 2, neg.VC = TRUE)                        

edfPlot(sim2)
VCPlot(sim2)

sim3 = simDataChisq(design =  design.df, 
            blk.str2 = "Run", blk.str1 = "Ani", 
            trt.str = "Tag + Trt", 
            gamma.R = c(0, 0.25, 1, 4, 100),                    
            gamma.A =c(10^((-8:8)/2)), 
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 2, neg.VC = FALSE, average = FALSE)                        

edfPlot(sim3)
VCPlot(sim3)

sim4 = simDataChisq(design =  design.df, 
            blk.str2 = "Run", blk.str1 = "Ani", 
            trt.str = "Tag + Trt", 
            gamma.R = c(0, 0.25, 1, 4, 100),                    
            gamma.A =c(10^((-8:8)/2)), 
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = FALSE)                        

pdf(width = 12, file = "CRD232Tag4Unadjusted.pdf")
edfPlot(sim4)
dev.off()

VCPlot(sim4)

sim = simData(design =  design.df, 
            blk.str2 = "Run", blk.str1 = "Ani", 
            trt.str = "Tag + Trt", 
            gamma.R = c(0, 0.1, 0.25, 1, 4, 100),                    
            gamma.A =c(10^((-8:8)/2)), 
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 1000, seed = 1,  neg.VC = TRUE)                        

sim1Run0 = simData(design =  design1.df, 
            blk.str2 = "Run", blk.str1 = "Ani", 
            trt.str = "Tag + Trt", 
            gamma.R = c(0), #0.1, 0.25, 1, 4, 100),                    
            gamma.A =c(10^((-8:8)/2)), 
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, neg.VC = FALSE)                        

plotEDF(sim1)

######################################################################################

old = proc.time()
design.df = optCRD(nTrt  = 2, bRep  = 3, tRep  = 2, nPlot = 4, iter  = 1000)
proc.time() - old

n = nrow(design.df)
nBlk = nlevels(design.df$Run)
nPlot = nlevels(design.df$Tag)
Run.mat = with(design.df, as.matrix(table(1:n, Run)))
Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))
Pb = projMat(Run.mat)
Pb1 = projMat(Tag.mat)

(ave.eff = test.CRD(X.trt = Trt.mat, (diag(n) - Pb) %*% (diag(n) - Pb1), Rep = n/ncol(Trt.mat))$ave.eff )

design1.df = optRBD(nTrt = 2, bRep  = 3, nCag  = 3, tRep  = 2, nPlot = 4, iter  = 3000, 
upperValue = ave.eff, confoundCag = FALSE)

design.summary.RBD(design1.df)

sim = simData1(design =  design1.df, 
            blk.str2 = "Run", blk.str1 = "Cag/Ani", 
            trt.str = "Tag + Trt", 
            gamma.R = c(0, 0.25, 1, 4, 100), 
            gamma.C = c(0, 0.25, 1, 4, 100),                    
            gamma.A =c(10^((-8:8)/2)), 
            VC.resid = 1,
            nS = 5, nVc = 4, nSim = 10000, neg.VC = FALSE)                        


str(sim$tempV)



plotEDF(sim, "Within RunBetweenCag:AniResidual") 

plotEDF1(sim,  "Within RunBetweenCag:AniResidual") 

  