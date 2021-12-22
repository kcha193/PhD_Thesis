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



design = data.frame(B = factor(rep(1:3, each = 10)),
                    R = factor(rep(1:10, 3)),
                    T = factor(rep(1:10, 3)),
                    V = factor(c( 1,2,1,3,4,3,2,5,5,4,
                                  4,5,2,4,5,1,3,3,1,2,
                                  5,1,3,5,2,4,4,2,3,1)))
                    
summaryAovOnePhase(design, blk.str = "R", trt.str = "V")

summaryAovTwoPhase(design,  blk.str2 = "B*T", blk.str1 = "R",
trt.str = "V")

                    
design = data.frame(Altitudes = factor(rep(1:3, each = 12)),
                    Soils = factor(rep(1:4, 9)),
                    Plants = factor(rep(1:12, 3)),
                    Benches = factor(rep(1:3, each = 4, time = 3)),
                    Viruse = factor(c(0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
                                      2, 2, 1, 1, 0, 0, 2, 2, 1, 1, 0, 0,
                                      1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 1, 1)))
design

summaryAovOnePhase(design, blk.str = "Soils/(Benches + Altitudes)", 
                            trt.str = "Soils * Viruse")



summaryAovTwoPhase(design, blk.str2 = "Benches/Plants*Altitudes", 
          blk.str1 = "Soils/(Benches + Altitudes)", trt.str = "Soils * Viruse")


#unblanced experiments


design = data.frame(B = factor(rep(1:4, each = 4)),
                    R = factor(c( 1,2,1,2,
                                  2,1,2,1,
                                  3,4,3,4,
                                  4,3,4,3)),
                    V = factor(c( 1,1,2,2,
                                  2,2,1,1,
                                  1,1,2,2,
                                  2,2,1,1)))

summaryAovOnePhase(design, blk.str = "R", trt.str = "V")

summaryAovTwoPhase(design,  blk.str2 = "B", blk.str1 = "R",
trt.str = "V")

summaryAovOnePhase(design, blk.str = "B", trt.str = "R")

summaryAovTwoPhase(design,  blk.str2 = "B", blk.str1 = "R",
trt.str = "V")


summaryAovTwoPhase(design[-1,],  blk.str2 = "B", blk.str1 = "R",
trt.str = "V")

summaryAovTwoPhase(design[-c(1:2),],  blk.str2 = "B", blk.str1 = "R",
trt.str = "V")

summaryAovOnePhase(design[-c(1:2),],  blk.str = "B", trt.str = "R")


y = rnorm(12)

summary(aov(y~ R + Error(B), design))

summary(aov(y~ V + Error(B + R ), design))

summaryAovTwoPhase(design,  blk.str2 = "B", blk.str1 = "R",
trt.str = "V", response = y)



#Biodiversity experiment


design = data.frame(Run = factor(rep(1:8, each = 8)), 
                    Method = factor(rep(1:2, each = 32)),
                    Sample = factor(rep(1:2, each = 32)),
                    Tillage = factor(rep(c( 1,1,1,1,
                                        1,1,1,1,
                                        2,2,2,2,
                                        2,2,2,2,
                                        2,1,2,1,
                                        2,1,2,1,
                                        1,2,1,2,
                                        1,2,1,2),2)), 
                    Block = factor(rep(c( 1,3,2,4,
                                          1,3,2,4,
                                          1,3,2,4,
                                          1,3,2,4,
                                          2,4,1,3,
                                          2,4,1,3,
                                          2,4,1,3,
                                          2,4,1,3),2)),
                    Plot = factor(rep(1:2, time = 4, each = 8)),
                    Depth = factor(rep(1:2, time = 8, each = 4)), 
                    Occasion = factor(rep(1:2, time = 16, each = 2)),
                    Interval = factor(rep(c(1,2,1,2,3,4,3,4), time = 8)),
                    Factions = factor(rep(1:2, time = 32)), 
                    Cluster = factor(rep(1:4, each = 16)),
                    Analysis = factor(rep(c(1,3,1,3,5,7,5,7,2,4,2,4,6,8,6,8), time = 4)))
                    
                    
summaryAovOnePhase(design,  blk.str = "Block/Plot/Sample + Depth", 
                                        trt.str = "(Method * Tillage) : Depth")


summaryAovOnePhase(design,  blk.str = "Occasion/Cluster", 
trt.str = "Block")


summaryAovTwoPhase(design,  blk.str2 = "Occasion/Cluster/Analysis", 
blk.str1 = "(Block/Plot/Sample + Depth)/Factions",
trt.str = "Method + Tillage")


summaryAovTwoPhase(design,  blk.str2 = "(Occasion/Interval) * Run", 
blk.str1 = "(Block/Plot/Sample + Depth)/Factions",
trt.str = "Method + Tillage")
