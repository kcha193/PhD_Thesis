#McIntyre's Design


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

require(infoDecompuTE)

#Phase 1 experiment
design1 <- local({ 
  Set = as.factor(rep(1:2, each = 16))
  TPla =  as.factor(rep(1:4,  each = 4, time =2) + rep(0:1, each = 16) * 4)
  TLfPos = as.factor(rep(1:4,  time = 8)) 
  Trt = as.factor(letters[c(1,2,3,4,
                            2,1,4,3,
                            3,4,1,2,
                            4,3,2,1,
                            1,3,4,2,
                            2,4,3,1,
                            3,1,2,4,
                            4,2,1,3)])
  Sam = as.factor(1:32)
  data.frame(Set, TPla, TLfPos, Trt, Sam)
})

design1


summaryAovOnePhase(design1, blk.str = "Set/(TPla*TLfPos)", trt.str = "Trt")

#Phase 2 experiment

design2 <- local({ 
  Squ = as.factor(rep(1:4, each = 32))
  APla = as.factor(rep(1:4,  time = 16, each = 2) + rep(0:3, each = 32) * 4)
  #APla = as.factor(rep(1:4,  time = 16, each = 2))

  ALfPos = as.factor(rep(1:4,  each = 8, time =4))
  Sam = c(1,17, 2,20, 3,18, 4,19,
          2,18, 1,19, 4,17, 3,20,
          3,19, 4,18, 1,20, 2,17,
          4,20, 3,17, 2,19, 1,18,
          
          5,23, 6,22, 7,24, 8,21,
          8,22, 7,23, 6,21, 5,24,
          7,21, 8,24, 5,22, 6,23,
          6,24, 5,21, 8,23, 7,22,
          
          9,28, 10,25, 11,27, 12,26, 
          10,27, 9,26, 12,28, 11,25,
          11,26, 12,27, 9,25, 10,28,
          12,25, 11,28, 10,26, 9,27,
          
          13,30, 14,31, 15,29, 16,32,
          16,31, 15,30, 14,32, 13,29,
          15,32, 16,29, 13,31, 14,30,
          14,29, 13,32, 16,30, 15,31)
  
  Set = design1$Set[match(Sam, design1$Sam)]
  TPla =design1$TPla[match(Sam, design1$Sam)]
  
  TLfPos = design1$TLfPos[match(Sam, design1$Sam)]
  Trt = design1$Trt[match(Sam, design1$Sam)]
  
  data.frame(Squ, APla, ALfPos, Set, TPla, TLfPos, Trt, Sam)
})

design2

summaryAovTwoPhase(design2, blk.str2 = "Squ/(APla*ALfPos)", 
                        blk.str1 = "Set/(TPla*TLfPos)", trt.str = "Trt")

summaryAovTwoPhase(design2, blk.str2 = "Squ/(APla*ALfPos)", 
                        blk.str1 = "Set*TPla*TLfPos", trt.str = "Trt")

summaryAovTwoPhase(design2, blk.str2 = "Squ+(APla*ALfPos)", 
                        blk.str1 = "Set + (TPla*(TLfPos/Set))", trt.str = "Trt")
                        
summaryAovTwoPhase(design2, blk.str2 = "Squ/(APla+ALfPos)", 
                        blk.str1 = "Set/(TPla*TLfPos)", trt.str = "Trt")
                        
################################################################################
#Efficiency factor

N = with(design2, table(interaction(Squ, APla, ALfPos), Trt))
A = 32*diag(4) - t(N) %*% N/2
C = diag(4) - matrix(1/4, nrow = 4, ncol = 4) 

sum(diag(C %*% t(C))) / (sum(diag(C %*% ginv(A) %*% t(C)))*32)

################################################################################

getVCs.onePhase(design2, trt.str = "TLfPos", 
blk.str = "Squ/(APla*ALfPos)",
var.comp = c("Squ:APla:ALfPos", "Squ:ALfPos", "Squ:APla"))


T1 = with(design2, sub("[12]", 1, TLfPos) )
T1 = as.numeric(sub("[34]", -1, T1))

T2 = with(design2, sub("[13]", 1, TLfPos) )
T2 = as.numeric(sub("[24]", -1, T2))

T3 = with(design2, sub("[14]", 1, TLfPos) )
T3 = as.numeric(sub("[23]", -1, T3))

getVCs.onePhase(design2, trt.str = "TLfPos", 
trt.contr = list(TLfPos = list(T1 = T1, T2 = T2, T3 = T3)),
blk.str = "Squ/(APla*ALfPos)",
var.comp = c("Squ:APla:ALfPos", "Squ:ALfPos", "Squ:APla"))



getVCs.onePhase(design2, trt.str = "Set/(TPla*TLfPos)", 
blk.str = "Squ/(APla*ALfPos)",
var.comp = c("Squ:APla:ALfPos", "Squ:ALfPos", "Squ:APla"))

getVCs.onePhase(design2, trt.str = "Trt", 
blk.str = "Squ/(APla*ALfPos)",
var.comp = c("Squ:APla:ALfPos", "Squ:ALfPos", "Squ:APla"))


summary.aov.twoPhase(design2, blk.str1 = "Set/(TPla*TLfPos)", 
blk.str2 = "Squ/(APla*ALfPos)",
trt.str = "Trt",
var.comp = c("Set:TPla:TLfPos", "Set:TLfPos", "Set:TPla",
"Squ:APla:ALfPos", "Squ:ALfPos", "Squ:APla"), latex = TRUE)  

##########################################################################################
d1.design = local({
   blk = factor(rep(1:12, each = 2))
   trt = factor(c(2,4,1,3,4,2,3,1,
                  3,2,4,1,1,4,2,3,
                  4,3,3,4,2,1,1,2))

   data.frame(blk,trt)
})

getVCs.onePhase(d1.design , trt.str = "trt", 
blk.str = "blk")

N = with(d1.design, table(blk, trt))

A = 6*diag(4) - t(N) %*% N/2

C = diag(4) - matrix(1/4, nrow = 4, ncol = 4) 

sum(diag(C %*% t(C))) / (sum(diag(C %*% ginv(A) %*% t(C)))*6)

d.design = local({
   blk = factor(rep(1:16, each = 2))
   trt = factor(c(1,1,2,2,3,3,4,4,
                  2,4,1,3,4,2,3,1,
                  3,2,4,1,1,4,2,3,
                  4,3,3,4,2,1,1,2))

   data.frame(blk,trt)
})

getVCs.onePhase(d.design , trt.str = "trt", 
blk.str = "blk")

N = with(d.design, table(blk, trt))

A = 8*diag(4) - t(N) %*% N/2

C = diag(4) - matrix(1/4, nrow = 4, ncol = 4) 

sum(diag(C %*% t(C))) / (sum(diag(C %*% ginv(A) %*% t(C)))*8)
