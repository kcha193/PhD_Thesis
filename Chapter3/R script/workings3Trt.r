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



design.df = optCRD(nTrt = 3, bRep  = 4, tRep  = 2, nPlot = 4, iter  = 1000)

nTrt = 3; nCag= 4; bRep  = 4; tRep  = 2; nPlot = 4;


design1.df = design.df

design.summary.RBD(design1.df)

Cage and Animal design:
  [,1] [,2] [,3] [,4]
[1,] "1A" "1B" "3A" "3B"
[2,] "1B" "1A" "3B" "3A"
[3,] "1C" "2A" "3C" "4A"
[4,] "2A" "1C" "4A" "3C"
[5,] "2B" "2C" "4B" "4C"
[6,] "2C" "2B" "4C" "4B"

design2.df = design.df

design.summary.RBD(design2.df)

[,1] [,2] [,3] [,4]
[1,] "1A" "1B" "1C" "2A"
[2,] "1B" "1A" "2A" "1C"
[3,] "2B" "2C" "3A" "3B"
[4,] "2C" "2B" "3B" "3A"
[5,] "3C" "4A" "4B" "4C"
[6,] "4A" "3C" "4C" "4B"













