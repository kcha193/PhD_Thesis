

install.packages(c("Rcpp", "inline", "RcppArmadillo"))  
 
library(MASS)
library(inline)
library(compiler)
library(Rcpp)
library(RcppArmadillo)


sourceDir <- function(path, trace = TRUE, ...) {
 
     for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        if (trace)
            cat(nm, ":")
        source(file.path(path, nm), ...)
        if (trace)
            cat("\n")
    }
}

sourceDir(path = "/home/kcha193/infoDecompuTE/R")
sourceDir(path = "/home/kcha193/optimTE/R")
load("/home/kcha193/powerAnalysis6trt.Rdata")

design.summary.RBD(design.df)

n = nrow(design.df)
nBlk = nlevels(design.df$Run)
nPlot = nlevels(design.df$Tag)
nCag = nlevels(design.df$Cag)

nAni = nlevels(design.df$Ani)
Run.mat = with(design.df, as.matrix(table(1:n, Run)))
Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
Cag.mat = with(design.df, as.matrix(table(1:n, Cag)))
Ani.mat = with(design.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))

fileName = paste("BIBDdesign_", ncol(Trt.mat), "trt_",  ncol(Ani.mat), "ani_", ncol(Cag.mat),
      "blk_",  ncol(Tag.mat), "tag_", ncol(Run.mat), "run.Rdata", sep = "")
                                                                               
fileName

save(design.df, file = fileName)
 ###################################################################################################
design.df = optBIBD(nTrt = 7, bRep  = 5, nCag  = 7, tRep  = 2, nPlot = 4, iter  = 100000, confoundCag = FALSE, resDF = 11)
design.df = optBIBD(nTrt = 6, bRep  = 4, nCag  = 6, tRep  = 2, nPlot = 4, iter  = 100000, confoundCag = FALSE, resDF = 11)
design.df = optBIBD(nTrt = 6, bRep  = 3, nCag  = 6, tRep  = 2, nPlot = 4, iter  = 100000, confoundCag = FALSE, resDF = 11)
design.df = optBIBD(nTrt = 6, bRep  = 2, nCag  = 6, tRep  = 2, nPlot = 4, iter  = 100000, confoundCag = FALSE, resDF = 11)



design.df = optBIBD(nTrt = 6, bRep  = 5, nCag  = 6, tRep  = 2, nPlot = 4, iter  = 100000, confoundCag = FALSE, resDF = 11)

design.df = optBIBD(nTrt = 6, bRep  = 6, nCag  = 6, tRep  = 2, nPlot = 4, iter  = 10000, confoundCag = FALSE, resDF = 11)

design.df = optBIBD(nTrt = 7, bRep  = 7, nCag  = 7, tRep  = 2, nPlot = 4, confoundCag = FALSE, iter = 10000)

design.df = optBIBD(nTrt =6, bRep  = 6, nCag  = 6, tRep  = 2, nPlot = 4, iter  = 10000, confoundCag = TRUE)

###################################################################################################


design.df = optRBD(nTrt = 4, bRep  = 6, nCag  = 6, tRep  = 2, nPlot = 4, iter  = 1000, confoundCag = FALSE)

design.df = optRBD(nTrt = 4, bRep  = 6, nCag  = 6, tRep  = 2, nPlot = 4, iter  = 1000, confoundCag = TRUE)

design.df = optRBD( nTrt = 8, bRep  = 9, nCag  = 9, tRep  = 2, nPlot = 4, iter = 5000, confoundCag = FALSE)

design.df = optRBD( nTrt = 8, bRep  = 10, nCag  = 10, tRep  = 2, nPlot = 4, iter = 10000, confoundCag = TRUE)

design.df = optRBD( nTrt = 8, bRep  = 10, nCag  = 10, tRep  = 2, nPlot = 4, iter = 10000, confoundCag = TRUE, upperValue = upper.ave.eff)


design.df = optRBD( nTrt = 4, bRep  = 9, nCag  = 3, tRep  = 2, nPlot = 4, iter  = 1000)

design.df = optRBD( nTrt = 4, bRep  = 9, nCag  = 3, tRep  = 2, nPlot = 4, iter  = 1000)

design.df = optRBD( nTrt = 3, bRep  = 8, nCag  = 4, tRep  = 2, nPlot = 4, iter  = 1000)

design.df = optRBD( nTrt = 6, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 4, iter  = 1000, confoundCag = TRUE)

###################################################################################################

design.df = optRBD( nTrt = 6, bRep  = 2, nCag  = 2, tRep  = 2, nPlot = 8, iter  = 3000, confoundCag = TRUE)

design.df = optRBD(nTrt = 6, bRep  = 6, nCag  = 6, tRep  = 2, nPlot = 8, iter  = 1000, confoundCag = FALSE)

design.df = optRBD( nTrt = 6, bRep  = 10, nCag  = 10, tRep  = 2, nPlot = 8, iter  = 1000, confoundCag = TRUE)

design.df = optRBD( nTrt = 4, bRep  = 5, nCag  = 5, tRep  = 2, nPlot = 8, iter  = 1000)

design.df = optRBD( nTrt = 6, bRep  = 8, nCag  = 8, tRep  = 2, nPlot = 8, iter  = 1000)
