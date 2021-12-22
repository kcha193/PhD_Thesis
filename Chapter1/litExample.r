




sourceDir <- function(path, trace = TRUE, ...) {
  library(MASS)
  library(inline)
  library(compiler) 
  library(formatR)
  library(Rcpp)
  library(RcppArmadillo)
  library(ggplot2)
  library(gridExtra)
  library(reshape)
  library(snowfall)
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
  }
}

sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R")

sourceDir(path = "C:/Users/Kevin/Dropbox/Packaging tools/packaging/infoDecompuTE/R")

 design = data.frame(Blk = factor(rep(1:4, each = 3)), 
                     Uit = factor(c(1,1,2,
                                    2,2,3,
                                    3,3,4,
                                    4,4,1)),
                     Trt = factor(c(1,1,2,
                                    2,2,1,
                                    1,1,2,
                                    2,2,1)), 
                     Sam = factor(1:12) )
design

summaryAovOnePhase(design.df = design, blk.str = "Uit", trt.str = "Trt")

summaryAovOnePhase(design.df = design, blk.str = "Blk", trt.str = "Uit")
    
summaryAovOnePhase(design.df = design, blk.str = "Blk", trt.str = "Trt")


summaryAovTwoPhase(design.df = design, blk.str2 = "Blk", blk.str1 = "Uit", 
    trt.str = "Trt")
    
summaryAovTwoPhase(design.df = design, blk.str2 = "Blk", blk.str1 = "Uit", 
    trt.str = "Sam")


contr1 =  c(1,1,-1,-1)[design$Uit]
contr2 =  c(1,-1,1,-1)[design$Uit]     #<- confounded with the treatment contrast
contr3 =   c(1,-1,-1,1)[design$Uit]


summaryAovOnePhase(design.df = design, blk.str = "Blk", trt.str = "Uit", 
trt.contr = list(Unit = list(contr1 = contr1, contr2 = contr2, contr3 = contr3)))



design = data.frame(Blk = factor(rep(1:4, each = 3)), 
                    Uit = factor(c(1,2,3,
                                    2,3,4,
                                    3,4,1,
                                    4,1,2)),
                    Trt = factor(c(1,2,1,
                                    2,1,2,
                                    1,2,1,
                                    2,1,2)), 
                    Sam = factor(1:12) )
design

summaryAovOnePhase(design.df = design, blk.str = "Uit", trt.str = "Trt")

summaryAovOnePhase(design.df = design, blk.str = "Blk", trt.str = "Uit")
    
summaryAovOnePhase(design.df = design, blk.str = "Blk", trt.str = "Trt")


summaryAovTwoPhase(design.df = design, blk.str2 = "Blk", blk.str1 = "Uit", 
    trt.str = "Trt")

n = nrow(design)
nBlk = nlevels(design$Blk)
nUit = nlevels(design$Uit)
nTrt = nlevels(design$Trt)

Blk.mat = with(design, as.matrix(table(1:n, Blk)))
Uit.mat = with(design, as.matrix(table(1:n, Uit)))
Trt.mat = with(design, as.matrix(table(1:n, Trt)))


(projMat(Trt.mat)- K(n)) %*% (identityMat(n) - projMat(Blk.mat)) %*% (projMat(Trt.mat)- K(n))
e = eigen((projMat(Trt.mat)- K(n)) %*% (identityMat(n) - projMat(Blk.mat) ) %*% (projMat(Trt.mat)- K(n)))

with(e, values[1] * vectors[,1] %*% t(vectors[,1] ))

(projMat(Trt.mat) - K(n)) * 8/9


(projMat(Blk.mat)- K(n)) %*% (projMat(Trt.mat)- K(n)) %*% (projMat(Blk.mat)- K(n))


(projMat(Blk.mat)- K(n)) %*% (projMat(Trt.mat)- K(n)) %*% (projMat(Blk.mat)- K(n)) * 9

(identityMat(n) - projMat(Blk.mat)) %*% (projMat(Trt.mat)- K(n)) %*% (identityMat(n) - projMat(Blk.mat)) * 9/8



tr((identityMat(n) - projMat(Blk.mat)) %*% (projMat(Uit.mat)- K(n)) %*% (identityMat(n) - projMat(Blk.mat)))/3


t(Uit.mat) %*% (identityMat(n) - projMat(Blk.mat)) %*% (projMat(Uit.mat)- K(n)) %*% (identityMat(n) - projMat(Blk.mat)) %*% Uit.mat

e = eigen((projMat(Trt.mat)- K(n)) %*% (identityMat(n) - projMat(Blk.mat) ) %*% (projMat(Trt.mat)- K(n)))

with(e, values[1] * vectors[,1] %*% t(vectors[,1] ))

(projMat(Trt.mat) - K(n)) * 8/9




(projMat(Uit.mat)- K(n)) %*% (identityMat(n) - projMat(Blk.mat) ) %*% (projMat(Uit.mat)- K(n))
e = eigen((projMat(Uit.mat)- K(n)) %*% (identityMat(n) - projMat(Blk.mat) ) %*% (projMat(Uit.mat)- K(n)))

with(e, vectors[,1:3] %*% t(vectors[,1:3] ))

with(e, values[1] * vectors[,1] %*% t(vectors[,1]) + 
        values[2] * vectors[,2] %*% t(vectors[,2]) +
        values[3] * vectors[,3] %*% t(vectors[,3]))


with(e, values[2] * vectors[,1] %*% t(vectors[,1] ))

(projMat(Trt.mat) - K(n)) * 8/9


(projMat(Trt.mat)- K(n)) %*% (projMat(Uit.mat)- K(n)) %*% (projMat(Trt.mat)- K(n))


(projMat(Trt.mat)- K(n))

N = t(Uit.mat) %*% Blk.mat
N %*% t(N)

Pb = projMat(Blk.mat)
blk.proj = (identityMat(n) - Pb)

fractions((eigen(t(Uit.mat) %*% (Pb - K(n)) %*% Uit.mat)$va)/(n/nUit))

fractions((eigen(t(Uit.mat) %*% blk.proj %*% Uit.mat)$va)/(n/nUit))
                                                                   

projMat(Uit.mat) %*% blk.proj  %*% projMat(Uit.mat) -  8/9 *projMat(Trt.mat)

projMat(Uit.mat)- projMat(Trt.mat)   


summaryAovOnePhase(design.df = design, blk.str = "Blk", trt.str = "Uit")
    


summaryAovTwoPhase(design.df = design, blk.str2 = "Blk", blk.str1 = "Uit", 
    trt.str = "Trt")
    

 design = data.frame(Blk = factor(rep(1:4, each = 3)), 
                     Uit = factor(c(1,2,3,
                                    4,1,2,
                                    3,4,1,
                                    2,3,4)),
                     Trt = factor(c(1,2,1,
                                    2,1,2,
                                    1,2,1,
                                    2,1,2)))

    
 design = data.frame(Blk = factor(rep(1:4, each = 3)), 
                     Uit = factor(c(1,2,3,
                                    2,3,4,
                                    3,4,1,
                                    4,1,2)),
                     Trt = factor(c(1,2,1,
                                    2,1,2,
                                    1,2,1,
                                    2,1,2)))

summaryAovOnePhase(design.df = design, blk.str = "Blk", trt.str = "Uit")
    


summaryAovTwoPhase(design.df = design, blk.str2 = "Blk", blk.str1 = "Uit", 
    trt.str = "Trt")
     