

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


#Athlete traing experiment
#9 training condition are to be investigated 
#from the combinations of 3 surface and 3 intensities

#12 athletes  divided into 3 lots for 3 months


design1 = data.frame( Months = factor(rep(1:4, each = 3*3)),
                      Athletes = factor(rep(1:3, times = 4, each = 3)),
                      Tests = factor(rep(1:3, times = 4*3)),
                      Condition = factor(rep(1:3, times = 4*3)))

summaryAovOnePhase(design.df = design1, blk.str = "Months/Athletes/Tests", 
                    trt.str = "Condition")



design2 = data.frame( Months = factor(rep(1:4, each = 3*3)),
                      Athletes = factor(rep(1:3, times = 4, each = 3)),
                      Tests = factor(rep(1:3, times = 4*3)),
                      Condition = factor(rep(1:3, times = 4*3)),
                      Surface = factor(rep(1:3, times = 4, each = 3)),
                      Batches = factor(rep(1:4, each = 3*3)),
                      Locations = factor(rep(1:9, times = 4)))


summaryAovTwoPhase(design.df = design2, 
                    blk.str2 = "Batches/Locations", 
                    blk.str1 = "Months /Athletes/Tests", 
                    trt.str = "Surface*Condition")


design3 = data.frame( Months = factor(rep(1:4, each = 3*3)),
                      Athletes = factor(rep(1:3, times = 4, each = 3)),
                      Tests = factor(rep(1:3, times = 4*3)),
                      Condition = factor(rep(1:3, times = 4*3)),
                      Surface = factor(rep(1:3, times = 4, each = 3)),
                      Batches = factor(rep(1:4, each = 3*3)),
                      Locations = factor(rep(1:3, times = 4, each = 3)),
                      Periods = factor(rep(1:3, times = 4*3)))

summaryAovOnePhase(design.df = design3, blk.str = "Batches/Periods/Locations", 
                    trt.str = "Surface*Condition")



summaryAovTwoPhase(design.df = design3, 
                    blk.str2 = "Batches/Periods/Locations", 
                    blk.str1 = "Months/Athletes/Tests", 
                    trt.str = "Surface*Condition")










