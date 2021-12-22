#workings
sourceDir <- function(path, trace = TRUE, ...) {
  library(MASS)
  library(inline)
  library(compiler) 
  library(formatR)
  library(ggplot2)
  library(gridExtra)
  library(reshape)
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
  }
}

sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R")
sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/optimTE/R")

sourceDir(path = "C:/Users/Kevin/Dropbox/Packaging tools/packaging/infoDecompuTE/R")
sourceDir(path = "C:/Users/Kevin/Dropbox/Packaging tools/packaging/optimTE/R")

fileName = paste("RBDdesign_", ncol(Trt.mat), "trt_",  ncol(Ani.mat), "ani_", ncol(Cag.mat),
      "blk_",  ncol(Tag.mat), "tag_", ncol(Run.mat), "run9.Rdata", sep = "")

fileName

save(design.df, file = fileName)

design.df = optCRD(nTrt = 2, bRep  = 4, tRep  = 2, nPlot = 4, iter  = 1000)


design.summary.CRD(design.df, FALSE)

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

design1.df = optRBD(nTrt = 2, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 4, 
			iter  = 1000, upperValue = ave.eff, confoundCag = FALSE)

design.summary.RBD(design1.df, FALSE)


design2.df = optRBD(nTrt = 2, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 4, 
iter  = 1000, upperValue = ave.eff, confoundCag = TRUE)


design.summary.RBD(design2.df)


design.df = optRBD(nTrt = 4, bRep  = 9, nCag  = 3, tRep  = 2, nPlot = 8, 
iter  = 3000, upperValue = ave.eff, confoundCag = FALSE)


design.df = optRBD(nTrt = 6, bRep  = 3, nCag  = 3, tRep  = 2, nPlot = 4, iter  = 10000, confoundCag = FALSE)

design.df = optRBD(nTrt = 8, bRep  = 10, nCag  = 10, tRep  = 2, nPlot = 8, 
            iter  = 2000, confoundCag = TRUE, upperValue = ave.eff)

summaryAovOnePhase(design1.df,  blk.str = "Run", trt.str = "Cag/Ani")

design.summary.RBD(design1.df)


design.summary.RBD(design2.df)

design.summary.RBD(design2.df)


#########################################################################################

old = proc.time()
design.df = optCRD(nTrt  = 6, bRep  = 4, tRep  = 2, nPlot = 4, iter  = 1000)
proc.time() - old    


design.summary.CRD(design.df)

design.summary.CRD(design.df, FALSE)


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

design1.df = optRBD(nTrt = 6, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 4, iter = 1000, confoundCag = TRUE, upperValue = ave.eff)

design.summary.RBD(design1.df, FALSE)


design2.df = optRBD(nTrt = 6, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 4, iter = 1000, confoundCag = FALSE, upperValue = ave.eff)

design.summary.RBD(design2.df, FALSE)



###########################################################################################
#RBD examples 

design.df = optRBD( nTrt = 4, bRep  = 2, nCag  = 2, tRep  = 2, nPlot = 8, iter  = 1000, upperValue = 3/5)


design.df = optRBD( nTrt = 8, bRep  = 8, nCag  = 2, tRep  = 2, nPlot = 8, iter  = 1000)

design.summary.RBD(design.df, FALSE)

design.summary.RBD(design.df)

###########################################################################################

design.df = optRBD(nTrt = 4, bRep  = 8, nCag  = 4, tRep  = 2, nPlot = 4, iter  = 1000)
   
design.df = optRBD(nTrt = 2, bRep  = 4, nCag  = 2, tRep  = 2, nPlot = 8, iter  = 1000)
   

##################################################################################

design.df = optRBD( nTrt = 8, bRep  = 10, nCag  = 10, tRep  = 2, nPlot = 4, resDF = 51, iter  = 10000)

design.summary.RBD(design.df, FALSE)

design.summary.RBD(design.df)

#################################################################################

design.df = optCRD( nTrt = 8, bRep  = 4, tRep  = 2, nPlot = 4, iter  = 1000)


fileName = paste("RBDdesign_", 8, "trt_", 32, "ani_", 4, "tag_", 
16, "run.Rdata", sep = "")
 
fileName

save(design.df, file = fileName)


trt.contr = test.CRD(X.trt = Trt.mat, (mI(n) - Pb) %*% (mI(n) - Pb1), Rep = n/ncol(Trt.mat))$e.vec
trt.contr = test.RBD(X.trt = Trt.mat, PP1, Rep = n/ncol(Trt.mat))$e.vec
            
trt.contr1 = trt.contr[,1][as.numeric(design.df$Trt)]
trt.contr2 = trt.contr[,2][as.numeric(design.df$Trt)]
trt.contr3 = trt.contr[,3][as.numeric(design.df$Trt)]
trt.contr4 = trt.contr[,4][as.numeric(design.df$Trt)]
trt.contr5 = trt.contr[,5][as.numeric(design.df$Trt)]
trt.contr6 = trt.contr[,6][as.numeric(design.df$Trt)]
trt.contr7 = trt.contr[,7][as.numeric(design.df$Trt)]


summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt",trt.contr = list(Tag = NA,
                                       Trt = list(  "1" = trt.contr1,
                                                    "2" = trt.contr2,
                                                    "3" = trt.contr3,
                                                    "4" = trt.contr4,
                                                    "5" = trt.contr5,
                                                    "6" = trt.contr6,
                                                    "7" = trt.contr7)), latex = TRUE,
                                       fixed.name = c("\\gamma", "\\tau_1", "\\tau_2",
                                                      "\\tau_3", "\\tau_4", "\\tau_5",
                                                      "\\tau_6","\\tau_7") )

design.summary.RBD(design.df)

summaryAovOnePhase(design.df, blk.str = "Cag/Ani", trt.str = "Trt", latex = TRUE)

summaryAovOnePhase(design.df, blk.str = "Run", trt.str = " Cag/Ani")

summaryAovTwoPhase(design14,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")


summaryAovTwoPhase(design14,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt",trt.contr = list(Tag = NA,
                                       Trt = list(  "1" = trt.contr1,
                                                    "2" = trt.contr2,
                                                    "3" = trt.contr3,
                                                    "4" = trt.contr4,
                                                    "5" = trt.contr5,
                                                    "6" = trt.contr6,
                                                    "7" = trt.contr7)), latex = TRUE,
                                       fixed.name = c("\\gamma", "\\tau_1", "\\tau_2",
                                                      "\\tau_3", "\\tau_4", "\\tau_5",
                                                      "\\tau_6", "\\tau_7") )
                                                    
                                                    
                                                   
                                                    
                                                    
                                                    