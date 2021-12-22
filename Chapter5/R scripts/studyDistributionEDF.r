#10/04/2013 10:08:58 a.m.
#Study the distrbution of the EDF upon deciding whether the mean or the median 
#of EDF's from the 

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


library(MASS)
library(ggplot2)
library(gridExtra)

design = design.df


aov.table = summaryAovTwoPhase(design.df = design, blk.str1 =   "Ani",
                          blk.str2 = "Run", trt.str = "Tag + Trt", response = rnorm(nrow(design)))

summaryAovTwoPhase(design.df = design, blk.str1 =  "Ani",
                          blk.str2 = "Run", trt.str = "Tag + Trt", response = rnorm(nrow(design)))

nSim = 10000

seed = 1

gamma.run = 0
gamma.ani = 1
VC.resid = 1

real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.run * VC.resid))

set.seed(seed)
simMS = apply(aov.table$A, 1, 
          function(x) sum(fracToNum(x[2:4]) * real.VC) * rchisq(nSim, fracToNum(x[1]))/fracToNum(x[1]))
nS
          
 apply(simMS, 2, mean)
 
old.time = proc.time()

neg.VC = TRUE
nS  = 4
nVc = 3  

REML.EDF = numeric(nS); LC.EDF = numeric(nS); REAL.EDF  = numeric(nS);
REML.VC = numeric(nVc); LC.VC = numeric(nVc); REAL.VC  = numeric(nVc);
  
 cat("gamma.ani = ", gamma.ani, "\n")
 pb <- txtProgressBar(min = 0, max = nSim, style = 3)
 counter = 0
 #simulation
 while( counter < nSim){
  setTxtProgressBar(pb, counter+1)
  
  real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.run * VC.resid))
    
  
  #aov.table$ANOVA[,"MS"] =  MS[counter+1, ]
  
  tmp = try(getVcEDF(aov.table = aov.table, MS = simMS[counter+1, ], row.MS = NA, true.VC = real.VC, neg.VC = neg.VC), TRUE)

  #print(tmp)
  #tmp
  if(class(tmp) =="try-error"){
  #stop()
  print(counter)
    counter = counter + 1
    next
  }
  
  tmpS = tmp$S
  tmpV = tmp$V
  
  REML.EDF  =  rbind(REML.EDF,     tmpS[,"REML.EDF"] )
  LC.EDF    =  rbind(LC.EDF  ,     tmpS[,"LC.EDF"]   )
  REAL.EDF  =  rbind(REAL.EDF,     tmpS[,"TRUE.EDF"] )
  
  REML.VC   =  rbind(REML.VC,      tmpV[,"REML.var.comp"]  )
  LC.VC     =  rbind(LC.VC  ,      tmpV[,"LC.var.comp"]    )
  REAL.VC   =  rbind(REAL.VC,      tmpV[,"TRUE.var.comp"]  )
 
  if((!is.na(tmpV[1,"REML.var.comp"])) && (tmpV[1,"REML.var.comp"]<0)) break 
  
  counter = counter + 1
}
close(pb)
proc.time() - old.time 
 
sum(is.na(REML.EDF[-1,1]))

which(is.na(REML.EDF[,1]))[1]

   
summaryEDF = cbind(apply(REAL.EDF[-1,], 2, function(x) median(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) mean(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) mean(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) median(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) median(x[which(!is.infinite(x))], na.rm = TRUE)))

colnames(summaryEDF) = c("True.EDF", "REML.EDF.mean", "LC.EDF.mean", 
                                     "REML.EDF.median", "LC.EDF.median")
summaryEDF

summaryVC = cbind(apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) median(x, na.rm = TRUE)))

colnames(summaryVC) = c("True.EDF", "REML.VC.mean", "LC.VC.mean", 
                                     "REML.VC.median", "LC.VC.median")
summaryVC
 
mean(apply(REML.VC[-1,], 1, function(x) x[1] + 2*x[2]), na.rm = TRUE)
mean(apply(LC.VC[-1,], 1, function(x) x[1] + 2*x[2]), na.rm = TRUE)
mean((simMS[, 2] - simMS[, 3])/2)
mean((simMS[, 8] - simMS[, 11])/2)

plot(apply(REML.VC[-1,], 1, function(x) x[1] + 2*x[2]), apply(LC.VC[-1,], 1, function(x) x[1] + 2*x[2]))
plot(REML.VC[-1,2], LC.VC[-1,2])
plot(REML.VC[-1,1], LC.VC[-1,1])

############################################################################################################  

Tag.mat = with(design, as.matrix(table(1:nrow(design), Tag)))
Trt.mat = with(design, as.matrix(table(1:nrow(design), Trt)))
   
X = cbind(1, Trt.mat, Tag.mat)

save.y = matrix(0, nrow = nSim, ncol = nrow(design))

seed = 1
old.time = proc.time()
gamma.run = 0
gamma.ani = 1
VC.resid = 1

true.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.run * VC.resid))

neg.VC = TRUE
REML.EDF = numeric(nS); LC.EDF = numeric(nS); REAL.EDF  = numeric(nS);
REML.VC = numeric(nVc); LC.VC = numeric(nVc); REAL.VC  = numeric(nVc);

nSim = 10000

cat("gamma.ani = ", gamma.ani, "\n")
pb <- txtProgressBar(min = 0, max = nSim, style = 3)
counter = 0
#simulation
while( counter < nSim){
  setTxtProgressBar(pb, counter+1)
   #set.seed(seed + counter)
                  
  set.seed(seed + counter);   run.eff = rnorm(nlevels(design$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
  set.seed(seed + counter+1); ani.eff = rnorm(nlevels(design$Ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
  set.seed(seed + counter+2); trt.eff = runif(nlevels(design$Trt), 0, 2)
  set.seed(seed + counter+3); tag.eff = runif(nlevels(design$Tag), 0, 1)
  set.seed(seed + counter+4); res.eff = rnorm(nrow(design), mean = 0, sd = sqrt(VC.resid))
  gm = 10       
  
  y = gm + with(design, run.eff[Run] + ani.eff[Ani] +  tag.eff[Tag] + trt.eff[Trt]) + res.eff

                     
  aov.table = summaryAovTwoPhase(design.df = design, blk.str1 =  blk.str1,
                    blk.str2 = blk.str2, trt.str = trt.str, 
                    response = y)
    
                    
  tmp = try(getVcEDF(aov.table = aov.table, row.MS = NA, true.VC = true.VC, neg.VC = neg.VC), TRUE)
                
  #tmp
                
  #print(tmp)   
  #tmp          
  if(class(tmp) =="try-error"){
    #stop()
    print(counter)
    counter = counter + 1
    next        
  }             
                
  tmpS = tmp$S  
  tmpV = tmp$V  
                
  REML.EDF  =  rbind(REML.EDF,     tmpS[,"REML.EDF"] )
  LC.EDF    =  rbind(LC.EDF  ,     tmpS[,"LC.EDF"]   )
  REAL.EDF  =  rbind(REAL.EDF,     tmpS[,"TRUE.EDF"] )
                
  REML.VC   =  rbind(REML.VC,      tmpV[,"REML.var.comp"]  )
  LC.VC     =  rbind(LC.VC  ,      tmpV[,"LC.var.comp"]    )
  REAL.VC   =  rbind(REAL.VC,      tmpV[,"TRUE.var.comp"]  )
                
  counter = counter + 1
  
  }        
  close(pb)
  proc.time() - old.time 

   dim(REML.EDF)
  
  sum(is.na(REML.VC[-1,1])     )
  
summaryEDF =cbind(apply(REAL.EDF[-1,], 2, function(x) median(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) mean(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) mean(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) median(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) median(x[which(!is.infinite(x))], na.rm = TRUE)))

colnames(summaryEDF) = c("True.EDF", "REML.EDF.mean", "LC.EDF.mean", 
                                     "REML.EDF.median", "LC.EDF.median")
summaryEDF

summaryVC = cbind(apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) median(x, na.rm = TRUE)))

colnames(summaryVC) = c("True.EDF", "REML.VC.mean", "LC.VC.mean", 
                                     "REML.VC.median", "LC.VC.median")
summaryVC
  
  
plot1 = qplot(REML.EDF[-1,1], geom="histogram", binwidth=0.01, xlab = "REML.EDF.BwRunBwAni")
plot2 = qplot(REML.EDF[-1,2], geom="histogram", binwidth=0.01, xlab = "REML.EDF.BwRunWiAni")
plot3 = qplot(REML.EDF[-1,3], geom="histogram", binwidth=0.1, xlab = "REML.EDF.WiRunBtwAni")
plot4 = qplot(REML.EDF[-1,4], geom="histogram", binwidth=0.01, xlab = "REML.EDF.WiRunWiAni")
 
plot5 = qplot(LC.EDF[-1,1], geom="histogram", binwidth=0.01, xlab = "LC.EDF.BwRunBwAni")
plot6 = qplot(LC.EDF[-1,2], geom="histogram", binwidth=0.01, xlab = "LC.EDF.BwRunWiAni")
plot7 = qplot(LC.EDF[-1,3], geom="histogram", binwidth=0.1, xlab = "LC.EDF.WiRunBtwAni")
plot8 = qplot(LC.EDF[-1,4], geom="histogram", binwidth=0.01, xlab = "LC.EDF.WiRunWiAni")

grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8,ncol=4)

apply(REML.EDF[-1,],2,range)
apply(LC.EDF[-1,],2,range)

pdf("CompareEDF.pdf")
grid.arrange( plot3 + geom_vline(xintercept = mean(REAL.EDF[-1,3]), col = "Red"), 
              plot7+ geom_vline(xintercept = mean(REAL.EDF[-1,3]), col = "Red"), ncol=2)
dev.off()

VC = data.frame(Method = factor(rep(c("LC", "REML"), each = nrow(REML.EDF[-1,]))),
                rbind(LC.VC[-1,], REML.VC[-1,]))


VC <- melt(VC)

 pg <- ggplot(VC, aes(value)) + 
 stat_density(geom = "histogram", position = "identity", 
 aes(colour = factor(Method))) + 
 facet_wrap(~ variable, scale = "fixed")  

  

 pg1 <- ggplot(VC, aes(value)) + 
 stat_density(geom = "histogram", position = "stack", 
 aes(colour = factor(Method))) + 
 facet_wrap(~ variable, scale = "fixed")  

  
 
grid.arrange(pg, pg1,  ncol=1)
 
 
plot1 = qplot(REML.VC[-1,1], geom="histogram", binwidth= VC.resid/10, xlab = "REML.VC.e")
plot2 = qplot(REML.VC[-1,2], geom="histogram", binwidth= VC.resid/10, xlab = "REML.VC.Ani")
plot3 = qplot(REML.VC[-1,3], geom="histogram", binwidth= VC.resid/2, xlab = "REML.VC.Run")
 
plot5 = qplot(LC.VC[-1,1], geom="histogram", binwidth= VC.resid/10, xlab = "LC.VC.e")
plot6 = qplot(LC.VC[-1,2], geom="histogram", binwidth= VC.resid/10, xlab = "LC.VC.Ani")
plot7 = qplot(LC.VC[-1,3], geom="histogram", binwidth= VC.resid/2, xlab = "LC.VC.Run")
 
grid.arrange(plot1, plot2, plot3, 
            plot5, plot6, plot7, ncol=3)
   
  
  
  
  
x1 = REML.VC[-1,1]
x1 = x1[which(!is.na(x1))]
x1 = cumsum(x1)/(1:length(x1))
plot1 = qplot(1:length(x1), x1, geom="path",  ylab = "REML.VC.e.mean", xlab = "Iteration")

x2 = REML.VC[-1,2]
x2 = x2[which(!is.na(x2))]
x2 = cumsum(x2)/(1:length(x2))
plot2 = qplot(1:length(x2), x2, geom="path",  ylab = "REML.VC.Ani.mean", xlab = "Iteration")

x3 = REML.VC[-1,3]
x3 = x3[which(!is.na(x3))]
x3 = cumsum(x3)/(1:length(x3))
plot3 = qplot(1:length(x3), x3, geom="path",  ylab = "REML.VC.Run.mean", xlab = "Iteration")

x4 = LC.VC[-1,1]
x4 = x4[which(!is.na(x4))]
x4 = cumsum(x4)/(1:length(x4))
plot4 = qplot(1:length(x4), x4, geom="path",  ylab = "LC.VC.e.mean", xlab = "Iteration")

x5 = LC.VC[-1,2]
x5 = x5[which(!is.na(x5))]
x5 = cumsum(x5)/(1:length(x5))
plot5 = qplot(1:length(x5), x5, geom="path",  ylab = "LC.VC.Ani.mean", xlab = "Iteration")

x6 = LC.VC[-1,3]
x6 = x6[which(!is.na(x6))]
x6 = cumsum(x6)/(1:length(x6))
plot6 = qplot(1:length(x6), x6, geom="path",  ylab = "LC.VC.Run.mean", xlab = "Iteration")

grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol=3)



x3 = REML.EDF[-1,3]
x3 = x3[which(!is.na(x3))]
x3 = cumsum(x3)/(1:length(x3))
plot3 = qplot(1:length(x3), x3, geom="path",  ylab = "REML.VC.Run.mean", xlab = "Iteration")


x6 = LC.EDF[-1,3]
x6 = x6[which(!is.na(x6))]
x6 = cumsum(x6)/(1:length(x6))
plot6 = qplot(1:length(x6), x6, geom="path",  ylab = "LC.VC.Run.mean", xlab = "Iteration")

grid.arrange( plot3, plot6, ncol=2)

x = REML.EDF[-1,2]
x = x[which(!is.na(x))]
x = cumsum(x)/(1:length(x))
qplot(1:length(x), x, geom="path")
abline(h= 2.028256)


x3 = sort(REML.VC[-1,2])
 mean( x3[1:(2*length(x3) - nSim)]) 


plot(1:200, x[1:200], type = "l")

x = REML.EDF[-1,2]
ind = which(!is.infinite(x))
x = LC.VC[-1,1]
x = x[ind]

x = cumsum(x)/(1:length(x))
plot(1:length(x), x, type = "l")
abline(h= 1)

x = REML.VC[-1,1]
x = cumsum(x)/(1:length(x))
plot(1:length(x), x, type = "l")


abline(h= 1)

x = sapply(lapply(1:length(x), function(y) 1:y), function(x) median(REML.VC[-1,3][x]))

plot(1:length(x), x, type = "l")
plot(1:100, x[1:100], type = "l")




x = REML.EDF[-1,3]
x = x[which(!is.na(x))]

x[which(x == max(x))]

x = x[-which(x == max(x))]
x[which(x == max(x))]

x = cumsum(x)/(1:length(x))
plot(1:length(x), x, type = "l")


x = LC.EDF[-1,3]
x = x[which(!is.na(x))]
x = cumsum(x)/(1:length(x))
plot(1:length(x), x, type = "l")

x = REML.EDF[-1,3]
x = x[which(!is.na(x))]

x[which(x == max(x))]

x = x[-which(x == max(x))]
x[which(x == max(x))]

x = cumsum(x)/(1:length(x))
plot(1:length(x), x, type = "l")







  
  pdf("filename.pdf")
grid.arrange(plot1, plot2, plot3, plot4,
            plot5, plot6, plot7, plot8,ncol=4)

dev.off()

#####################################################################################

aov.table = summaryAovTwoPhase(design.df = design.df, blk.str1 =  "Cag/Ani",
                          blk.str2 = "Run", trt.str = "Tag + Trt", response = rnorm(nrow(design.df)))

nSim = 10000

gamma.run = 0
gamma.cag = 0.01
gamma.ani = 10
VC.resid = 1
 
real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid),(gamma.run * VC.resid))

simMS = apply(aov.table$A, 1, 
          function(x) sum(fracToNum(x[2:5]) * real.VC) * rchisq(nSim, fracToNum(x[1]))/fracToNum(x[1]))
          
 
old.time = proc.time()

neg.VC = TRUE

nS =  6; nVc = 4

REML.EDF = numeric(nS); LC.EDF = numeric(nS); REAL.EDF  = numeric(nS);
REML.VC = numeric(nVc); LC.VC = numeric(nVc); REAL.VC  = numeric(nVc);
  
 cat("gamma.ani = ", gamma.ani, "\n")
 pb <- txtProgressBar(min = 0, max = nSim, style = 3)
 counter = 0
 #simulation
 while( counter < nSim){
  setTxtProgressBar(pb, counter+1)
      
  
  #aov.table$ANOVA[,"MS"] =  MS[counter+1, ]
  
  tmp = try(getVcEDF(aov.table = aov.table, MS = simMS[counter+1, ], row.MS = NA, true.VC = real.VC, neg.VC = neg.VC), TRUE)
   # tmp = try(getVcEDF.modified(aov.table = aov.table, row.MS = NA, true.VC = real.VC, neg.VC = neg.VC), TRUE)

  #print(tmp)
  #tmp
  if(class(tmp) =="try-error"){
  #stop()
  #print(counter)
    counter = counter + 1
    next
  }
  
  tmpS = tmp$S
  tmpV = tmp$V
  
  REML.EDF  =  rbind(REML.EDF,     tmpS[,"REML.EDF"] )
  LC.EDF    =  rbind(LC.EDF  ,     tmpS[,"LC.EDF"]   )
  REAL.EDF  =  rbind(REAL.EDF,     tmpS[,"TRUE.EDF"] )
  
  REML.VC   =  rbind(REML.VC,      tmpV[,"REML.var.comp"]  )
  LC.VC     =  rbind(LC.VC  ,      tmpV[,"LC.var.comp"]    )
  REAL.VC   =  rbind(REAL.VC,      tmpV[,"TRUE.var.comp"]  )
  
  counter = counter + 1
}
close(pb)
proc.time() - old.time 
 

       
summaryEDF = cbind(apply(REAL.EDF[-1,], 2, function(x) median(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) median(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) median(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) mean(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) mean(x[which(!is.infinite(x))], na.rm = TRUE)))

colnames(summaryEDF) = c("True.EDF", "REML.EDF.median", "LC.EDF.median", 
                                     "REML.EDF.mean", "LC.EDF.mean")
summaryEDF

summaryVC = cbind(apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) median(x, na.rm = TRUE)))

colnames(summaryVC) = c("True.EDF", "REML.VC.mean", "LC.VC.mean", 
                                     "REML.VC.var", "LC.VC.var")
summaryVC
       