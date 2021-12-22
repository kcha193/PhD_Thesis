

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

setwd("C:/Users/kcha193/My Dropbox/Chapter 6 - Estimating VC and EDF/Chris' design")

design = read.csv("Oats.csv")
str(design)
design$Blocks = as.factor(design$Blocks )
design$Wplots = as.factor(design$Wplots )
design$Subplots = as.factor(design$Subplots)

design$Variety = as.factor(as.numeric(design$Variety ))
design$Nitrogen = as.factor(as.numeric(design$Nitrogen ))

 
 aov.table = summaryAovOnePhase(design, blk.str = "Blocks/Wplots/Subplots", 
        trt.str = "Variety * Nitrogen", response = as.numeric(design$Y))
  
 
sqrt(2 * 601.33056  / 24)

sqrt(2 * 177.08333 / 18)

sqrt(2 * 177.08333 / 6)

MS2 = 601.3
MS3 = 177.1
 
sW = (MS2 - MS3)/4 
sE = MS3

sqrt(2 *(sW + sE)/ 6)

nu = (177.08333 + 106.0618)^2
den =  177.08333 ^2/(45) +  (601.33056/4)^2/(10) + (177.08333 /4)^2/(45)

 nu/den
   
 3/4 * MS3 +  MS2/4 
  
den =   (3/4 * MS3)^2 / 45 +   (MS2/4)^2 / 10 
 
  
  