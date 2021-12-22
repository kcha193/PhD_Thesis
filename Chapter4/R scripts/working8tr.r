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


design.summary.RBD(design.df)

sim = simData1(design =  design.df, 
            blk.str2 = "Run", blk.str1 = "Cag/Ani", 
            trt.str = "Tag + Trt", 
            gamma.R = c(0), #0.25 , 1, 4, 100), 
            gamma.C = c(0), #0.25, 1, 4, 100),                    
            gamma.A = c(10^((-10:10)/2)), 
            VC.resid = 1,
            nS = 4, nVc = 4, nSim = 100, neg.VC = TRUE)                        


str(sim$tempV)



plotEDF1(sim)

plotEDF2(sim)

plotEDF1(sim,  "Within RunBetweenCag:AniResidual") 
##################################################################################################
design.summary.CRD(designCRD, FALSE)

design.summary.RBD(design16, FALSE)

design.summary.RBD(design14, FALSE)


###############################################################################################
design.df$Ani =  with(design.df, interaction(Cag, Ani))
design14$ani =  with(design14, interaction(Cag, Ani))
design16$ani =  with(design16, interaction(Cag, Ani))

(anovaTable = summaryAovTwoPhase(design.df = design.df ,blk.str2 = "Run", 
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)
 
 (anovaTable14 = summaryAovTwoPhase(design.df = design14 ,blk.str2 = "Run", 
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)
 
 (anovaTable16 = summaryAovTwoPhase(design.df = design16 ,blk.str2 = "Run", 
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)
 
  (anovaTableCRD = summaryAovTwoPhase(design.df = designCRD ,blk.str2 = "Run", 
              blk.str1 = "Ani", trt.str = "Tag + Trt")$A)
                
                  
VC.resid = 1

gamma.run=  10
gamma.cag =  5
 gamma.ani =  2 #  c(10^((-4:4)/2))

nSim = 100000


# create progress bar
for (i in 1:length(gamma.R)) {
    gamma.run = gamma.R[i]
    
    cat("gamma.run = ", gamma.run, "\n")
    for (k in 1:length(gamma.C)) {
        gamma.cag = gamma.C[k]
        cat("gamma.cag = ", gamma.cag, "\n")
        
        for (j in 1:length(gamma.A)) {
            
            gamma.ani = gamma.A[j]
            
            
            MS14 = rep(0, nrow(anovaTable14))
            MS16 = rep(0, nrow(anovaTable16))
            #MSCRD = rep(0, nrow(anovaTableCRD))
            
            cat("gamma.ani = ", gamma.ani, "\n")
            pb <- txtProgressBar(min = 0, max = nSim, style = 3)
            counter = 0
            # simulation
            
            
            while (counter < nSim) {
                setTxtProgressBar(pb, counter + 1)
                
                run.eff = rnorm(nlevels(design.df$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
                cag.eff = rnorm(nlevels(design.df$Cag), mean = 0, sd = sqrt(gamma.cag * VC.resid))
                ani.eff = rnorm(nlevels(design.df$Ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
               trt.eff = c(1,1,1,1,4,4,4,4) #runif(nlevels(design.df$Trt), 0, 10)
                tag.eff = runif(nlevels(design.df$Tag), 0, 1)
                res.eff = rnorm(nrow(design.df), mean = 0, sd = sqrt(VC.resid))
                gm = 10
                
                
                # yCRD = gm + with(designCRD, run.eff[Run] + ani.eff[Ani] + 
                #        tag.eff[Tag] + trt.eff[Trt]) + res.eff
                        
                y14 = gm + with(design14, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
                        tag.eff[Tag] + trt.eff[Trt]) + res.eff
                
                y16 = gm + with(design16, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
                        tag.eff[Tag] + trt.eff[Trt]) + res.eff
                
               #  aov.tableCRD = summaryAovTwoPhase(design.df = designCRD, blk.str2 = "Run", 
               # blk.str1 = "Ani", trt.str = "Tag + Trt", 
               #   response = yCRD)
             
                aov.table14 = summaryAovTwoPhase(design.df = design14, blk.str2 = "Run", 
                blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
                  response = y14)
                
                aov.table16 = summaryAovTwoPhase(design.df = design16, blk.str2 = "Run", 
                blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
                  response = y16)
                
                #aov.table14
                
                #summary(aov(y14~ Tag + Trt+ Error(Run + Cag/Ani), data = design14))
   
                #aov.table16

                #summary(aov(y16~ Tag + Trt+ Error(Run + Cag/Ani), data = design16))
                   
                
                #summary.lm(aov(y14~ Tag + Trt+ Error(Run + Cag/Ani), data = design14)$"Cag:Ani")
                #summary.lm(aov(y16~ Tag + Trt+ Error(Run + Cag/Ani), data = design16)$"Cag:Ani")

             #   MSCRD = cbind(MSCRD, as.numeric(aov.tableCRD$A[, "MS"]))
                
                MS14 = cbind(MS14, as.numeric(aov.table14$A[, "MS"]))
                MS16 = cbind(MS16, as.numeric(aov.table16$A[, "MS"]))
                
                counter = counter + 1
            }
            close(pb)
            
          #  MSCRD = MSCRD[, -1]
        
            MS14 = MS14[, -1]
            
            MS16 = MS16[, -1]
                  
                  
                  
                     
            temp14 = rbind(temp14, cbind(gamma.run, gamma.cag, gamma.ani, mean(MS14[14, ]), mean(MS14[11, ] - MS14[14, 
                ])/2, mean(MS14[8, ] - MS14[11, ])/16, mean(MS14[4, ] - MS14[14, ])/4, mean(MS14[11, ])))
            
            temp16 = rbind(temp16, cbind(gamma.run, gamma.cag, gamma.ani, mean(MS16[13, ]), 
                    mean(MS16[10, ] - MS16[13,])/2, 
                    mean(MS16[2, ] - MS16[5, ] - MS16[10, ] + MS16[13, ])/16, 
                    mean(MS16[5, ] - MS16[13, ])/4, mean(MS16[10, ])))
            
        }
        
    }
} 


#######################################################################################
aov.table14 = summaryAovTwoPhase(design.df = design14, blk.str2 = "Run", blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
                  response = rnorm(nrow(design14)))

aov.table16 = summaryAovTwoPhase(design.df = design16, blk.str2 = "Run", blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
                  response = rnorm(nrow(design16)))

 gamma.run = 100
 gamma.cag = 10
 gamma.ani =5
 
 real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag *  VC.resid), (gamma.run * VC.resid))
  
 MS14.ChiSqu = t(suppressWarnings(apply(aov.table14$A, 1, 
          function(x) sum(fracToNum(x[2:5]) * real.VC) * rchisq(100000, fracToNum(x[1]))/fracToNum(x[1]))))

 MS16.ChiSqu = t(suppressWarnings(apply(aov.table16$A, 1, 
          function(x) sum(fracToNum(x[2:5]) * real.VC) * rchisq(100000, fracToNum(x[1]))/fracToNum(x[1]))))
 

library(ggplot2)
library(gridExtra)

plot1 = qplot(MS14.ChiSqu[11,]/8, geom="histogram", binwidth=.1, xlab = "Treatment variance with 14 DF")
plot2 = qplot(MS16.ChiSqu[10,]/8, geom="histogram", binwidth=.1, xlab = "Treatment variance with 16 DF")
plot3 = qplot((MS14.ChiSqu[11,]- MS14.ChiSqu[14,])/2, geom="histogram", binwidth=.1, xlab = "VC of animals within cages with 14 DF")
plot4 = qplot((MS16.ChiSqu[10,] -  MS16.ChiSqu[13,])/2, geom="histogram", binwidth=.1, xlab = "VC of animals within cages with 16 DF")
plot5 = qplot(MS14.ChiSqu[14,], geom="histogram", binwidth=.1, xlab = "Base VC with 14 DF")
plot6 = qplot(MS16.ChiSqu[13,], geom="histogram", binwidth=.1, xlab = "Base VC with 16 DF")

grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol=2)

mean(MS14.ChiSqu[11,]/8)
mean(MS16.ChiSqu[10,]/8)

mean((MS14.ChiSqu[11,]- MS14.ChiSqu[14,])/2)
mean((MS16.ChiSqu[10,] -  MS16.ChiSqu[13,])/2)

mean(MS14.ChiSqu[14,])
mean(MS16.ChiSqu[13,])

######################################################################################

library(ggplot2)
library(gridExtra)

plot1 = qplot(MS14[11,]/8, geom="histogram", binwidth=.1, xlab = "Treatment variance with 14 DF")
plot2 = qplot(MS16[10,]/8, geom="histogram", binwidth=.1, xlab = "Treatment variance with 16 DF")
plot3 = qplot((MS14[11,]- MS14[14,])/2, geom="histogram", binwidth=.5, xlab = "VC of animals within cages with 14 DF")
plot4 = qplot((MS16[10,] -  MS16[13,])/2, geom="histogram", binwidth=.5, xlab = "VC of animals within cages with 16 DF")
plot5 = qplot(MS14[14,], geom="histogram", binwidth=.1, xlab = "Base VC with 14 DF")
plot6 = qplot(MS16[13,], geom="histogram", binwidth=.1, xlab = "Base VC with 16 DF")

grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol=2)

mean(MS16[10,]/(8 * 0.84))
mean(MS14[11,]/(8 * 14259/16780))

mean((MS16[10,] -  MS16[13,])/2)
mean((MS14[11,]- MS14[14,])/2)

mean(MS16[9,]/MS16[10,])
mean(MS14[10,]/MS14[11,])

sum(1 - pf(MS16[9,]/MS16[10,], 7, 16) < 0.05)/ nSim
sum(1 - pf(MS14[10,]/MS14[11,], 7, 14) < 0.05)/ nSim

sum(p.adjust(1 - pf(MS16[9,]/MS16[10,], 7, 16), "fdr") < 0.05)/nSim
sum(p.adjust(1 - pf(MS14[10,]/MS14[11,], 7, 14), "fdr") < 0.05)/ nSim

 real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag *  VC.resid), (gamma.run * VC.resid))

VCDist = function(interval){
  
  plot1 = qplot(MS14[14,interval], geom="histogram", binwidth=.1, xlab = "VC of measurement error with 14 DF")
  plot2 = qplot(MS16[13,interval], geom="histogram", binwidth=.1, xlab = "VC of measurement error with 16 DF")

  
  plot4 = qplot((MS14[11,interval]- MS14[14,interval])/2, geom="histogram", binwidth=.1, xlab = "VC of animals within cages with 14 DF")
  plot5 = qplot((MS16[10,interval] -  MS16[13,interval])/2, geom="histogram", binwidth=.1, xlab = "VC of animals within cages with 16 DF")

  
  plot7 = qplot((MS14[8,interval] -MS14[11,interval])/16, geom="histogram", binwidth=.5, xlab = "VC of cages with 14 DF")
  plot8 = qplot((MS16[2,interval] - MS16[5,interval] - MS16[10,interval] +  MS16[13,interval])/16, 
                geom="histogram", binwidth=.5, xlab = "VC of cages with 16 DF")

  plot9 = qplot((MS14[4,interval] -  MS14[14,interval])/4, geom="histogram", binwidth=.5, xlab = "VC of runs with 14 DF")
  plot10 = qplot((MS16[5,interval] -  MS16[13,interval])/4, geom="histogram", binwidth=.5, xlab = "VC of runs with 16 DF")
  
  grid.arrange(plot2, plot1,  plot5,  plot4,plot8, plot7,  plot10,  plot9,ncol=2)
  
   
  result = data.frame("TRUE.value" =c(VC.resid, (gamma.ani * VC.resid), (gamma.cag *  VC.resid), (gamma.run * VC.resid)),
  DF16mean  =c(mean(MS16[13,interval]), 
              mean((MS16[10,interval] -  MS16[13,interval])/2), 
              mean((MS16[2,interval] - MS16[5,interval] - MS16[10,interval] +  MS16[13,interval])/16), 
              mean((MS16[5,interval] -  MS16[13,interval])/4)),
 DF14mean  = c(mean(MS14[14,interval]),  
              mean((MS14[11,interval]- MS14[14,interval])/2), 
              mean((MS14[8,interval] -MS14[11,interval])/16), 
               mean((MS14[4,interval] -  MS14[14,interval])/4)))
  
  
  
    rownames(result) = c("Measurement error", "Animals within cages", "Cages", "Runs")
  
  return(result)
}

VCDist(sample( 1:nSim, nSim))


getwd()
pdf(file = "VC8trt.pdf")
VCDist( sample(1:nSim, nSim))
dev.off()



trtDist = function(interval){

    plot1 = qplot(MS14[11,interval]/(8* 14259/16780 ), geom="histogram", binwidth=.1, xlab = "Treatment variance with 14 DF")
   plot2 = qplot(MS16[10,interval]/(8* 0.84 ), geom="histogram", binwidth=.1, xlab = "Treatment variance with 16 DF")
 
   
  plot7 = qplot(MS14[10,interval]/MS14[11,interval], geom="histogram", binwidth=.5, xlab = "F-ratio with 14 DF")
   plot8 = qplot(MS16[9,interval]/MS16[10,interval], geom="histogram", binwidth=.5, xlab = "F-ratio with 16 DF")
  
   
   plot4 = qplot(p.adjust(1 - pf(MS14[10,interval]/MS14[11,interval], 7, 14), "fdr"), geom="histogram", binwidth=.05, xlab = "Adjusted p-value with 14 DF")
   
  plot5 = qplot(p.adjust(1 - pf(MS16[9,interval]/MS16[10,interval], 7, 16), "fdr"), geom="histogram", binwidth=.05, xlab = "Adjusted p-value with 16 DF")   

  grid.arrange( plot2, plot1, plot8,  plot7, plot5, plot4,   ncol=2)
 

  result = data.frame(  
  DF16mean  =c(mean(MS16[10,interval]/(8* 0.84 )), 
          mean(MS16[9,interval]/MS16[10,interval]),
          sum(p.adjust(1 - pf(MS16[9,interval]/MS16[10,interval], 7, 16), "fdr") < 0.05)/length(interval)),
  DF14mean  = c(mean(MS14[11,interval]/(8* 14259/16780) ), 
          mean(MS14[10,interval]/MS14[11,interval]),
          sum(p.adjust(1 - pf(MS14[10,interval]/MS14[11,interval], 7, 14), "fdr") < 0.05)/ length(interval)))
 
  
  rownames(result) = c("Treatment variance", "F-ratio", "Power of test")
  
  return(result)
}


getwd()
pdf(file = "trt8trt.pdf")
trtDist(sample( 1:nSim, nSim))
dev.off()


###################################################################################################################

#Compare the SED one a single data set  

#Pairwise comparison
contr.e = contrCompare(design16, pairwiseContr(design16))

contrastMat16 = matrix(0, ncol= 8, nrow = 8)

for( i in 1:length(contr.e)){   
  contrastMat16[states[i,2], states[i,1]] = contr.e[i]
}

trtVar = rowMeans(sapply(MS16[10,], function(x) sqrt(2*x/(8* contr.e)))) 

for( i in 1:length(trtVar)){   
  contrastMat16[states[i,1], states[i,2]] = trtVar[i]
}

contrastMat16


contr.e = contrCompare(design14, pairwiseContr(design14))

contrastMat14 = matrix(0, ncol= 8, nrow = 8)

for( i in 1:length(contr.e)){   
contrastMat14[states[i,2], states[i,1]] = contr.e[i]
}

trtVar = rowMeans(sapply(MS14[11,], function(x) sqrt(2*x/(6* contr.e)))) 


for( i in 1:length(trtVar)){   
  contrastMat14[states[i,1], states[i,2]] = trtVar[i]
}

contrastMat14


sum(contrCompare(design16, pairwiseContr(design16)) > contrCompare(design14, pairwiseContr(design14)))
sum(contrCompare(design16, pairwiseContr(design16)) < contrCompare(design14, pairwiseContr(design14)))



   
contr1 = c(1, 0, 0, -1,0,0,0,0)
contr2 =c(0, 1, 0,0,0,0,-1,0)
contr3 = c( 0, 0,1, 0, 0, 0,0,-1)
contr4 =c( 0, 0,0, 0, 1, -1,0,0)
contr5 =  abs(contr1 ) - abs(contr2)
contr6 =  abs(contr3 ) - abs(contr4)
contr7 = abs(contr5 ) - abs(contr6)


contr = cbind(contr1, contr2, contr3, contr4, contr5, contr6, contr7)/2  
contr

contrCompare(design16, contr)
contrCompare(design14, contr)

sqrt(2*as.numeric(aov.table16[10,"MS"])/(8* contrCompare(design16, contr))  )
  
  contrCompare(design6, basicContr(design6)[,-6])
  
  
(aov.table14 = summaryAovTwoPhase(design.df = design14, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
      response = y14)$A)
   
 aov14.fit = aov(y14~  Tag + Trt + Error(Run + Cag/Ani), design14)
  summary(aov14.fit)

  summary.lm(aov14.fit$"Cag:Ani")
  
  

 











##################################################################################################################

plot1 = qplot(rf(10000, 7, 14), geom="histogram", binwidth=.1, xlab = "Random numbers from F-dist with 14 DF")
plot2 = qplot(rf(10000, 7, 15), geom="histogram", binwidth=.1, xlab = "Random numbers from F-dist with 16 DF")
grid.arrange(plot1, plot2, ncol=2)

sum((1- pf(MS16[9,]/MS16[10,], 7, 16)) < 0.05)
sum((1- pf(MS14[10,]/MS14[11,], 7, 14)) < 0.05)

############################################################################################
MS14 = MS14[,-1]

c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid), (gamma.run * VC.resid))


 mean(MS14[14,])
 mean(MS14[11,] -  MS14[14,])/2
mean(MS14[8,] -  MS14[11,])/16 
 mean(MS14[4,] -  MS14[14,])/4
  mean(MS14[11,])
  
  
#########################################################################################
> c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid), (gamma.run * VC.resid))
[1]   1   5  10 100
> 
> anovaTable
                   DF e Cag:Ani Cag Run
Between Run                            
   Between Cag:Ani                     
      Trt          7  1 2       0   4  
   Within Cag.Ani  8  1 0       0   4  
Within Run                             
   Between Cag                         
      Tag          1  1 2       16  0  
      Residual     2  1 2       16  0  
   Between Cag:Ani                     
      Trt          7  1 2       0   0  
      Residual     14 1 2       0   0  
   Within Cag.Ani                      
      Tag          2  1 0       0   0  
      Residual     22 1 0       0   0  
>   
> mean(MS[14,])
[1] 0.9964707
> 
> mean(MS[11,] -  MS[14,])/2
[1] 4.964581                                  
> 
> mean(MS[7,] -  MS[11,])/16
[1] 10.11867
> 
> mean(MS[4,] -  MS[14,])/4
[1] 100.3515
> 
> 
> mean(MS[11,])
[1] 10.92563


############################################################################################
MS16 = MS16[,-1]
c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid), (gamma.run * VC.resid))

mean(MS16[13,])
 mean(MS16[10,] -  MS16[13,])/2
 mean(MS16[2,] -  MS16[5,]-  MS16[10,] +   MS16[13,])/16
  mean(MS16[5,] -  MS16[13,])/4
mean(MS16[10,])
  
############################################################################################

> c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid), (gamma.run * VC.resid))
[1]   1   5  10 100
> 
> anovaTable
                   DF e Cag:Ani Cag Run
Between Run                            
   Between Cag     3  1 2       16  4  
   Between Cag:Ani                     
      Trt          4  1 2       0   4  
   Within Cag.Ani  8  1 0       0   4  
Within Run                             
   Between Cag:Ani                     
      Tag          1  1 2       0   0  
      Trt          7  1 2       0   0  
      Residual     16 1 2       0   0  
   Within Cag.Ani                      
      Tag          2  1 0       0   0  
      Residual     22 1 0       0   0  
>   
> mean(MS[13,])
[1] 1.001878
> 
> mean(MS[10,] -  MS[13,])/2
[1] 5.014447
> 
> mean(MS[2,] -  MS[5,]-  MS[10,] +   MS[13,])/16
[1] 9.813845
> 
> mean(MS[5,] -  MS[13,])/4
[1] 100.8395
> 
> mean(MS[10,])
[1] 11.03077



mean(MS[13])
 mean(MS[10] -  MS[13])/2
 mean(MS[2] -  MS[5]-  MS[10] +   MS[13])/16
  mean(MS[5] -  MS[13])/4
mean(MS[10])


  
 VC.resid = 1
gamma.run = 10
  gamma.cag = 0
   gamma.ani =  100
 
 

 run.eff = rnorm(nlevels(design.df$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
  cag.eff = rnorm(nlevels(design.df$Cag), mean = 0, sd = sqrt(gamma.cag  * VC.resid))
  ani.eff = rnorm(nlevels(design.df$ani), mean = 0, sd = sqrt( gamma.ani * VC.resid))
  trt.eff = runif(nlevels(design.df$Trt), 0, 2)
  tag.eff = runif(nlevels(design.df$Tag), 0, 1)
  res.eff = rnorm(nrow(design.df), mean = 0, sd = sqrt(VC.resid))
  gm = 10
  
  y =   gm + with(design.df, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] +  tag.eff[Tag] + trt.eff[Trt]) + res.eff
  
  aov.table = summaryAovTwoPhase(design.df = design.df ,blk.str2 = "Run", blk.str1 = "Cag/Ani", 
                 trt.str = "Tag + Trt", response = y)
    
  m1.asreml = asreml(y ~ Tag + Trt,
                           random = ~ Run + Cag/Ani,
                           data = design.df)
 
 summary(m1.asreml)$varcomp
 wald(m1.asreml, denDF="algebraic")
 
  getVcEDF(aov.table = aov.table, row.MS = NA, 
          true.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid), (gamma.run * VC.resid)), 
            neg.VC = FALSE)
 
 
#######################################################################################
basicContr = function(design.df){

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

    C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
    C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
    cage.Rep = n/nCag

    C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))

    Pb = projMat(Run.mat)
    Pb1 = projMat(Tag.mat)

    blk.proj = (identityMat(n) - Pb) %*% (identityMat(n) - Pb1)


     
      info.mat = matMulti(blk.proj, Ani.mat, C.cage)
  
      if (any(abs(info.mat) > 1e-07)) {
          PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)),
              C.cage, Ani.mat)
  
          PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
      } else {
          PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani,
              Ani.mat)
  
      }
  
      test.RBD(X.trt = Trt.mat, PP1, Rep = n/ncol(Trt.mat))$e.vec

}


##########################################################################################
contrCompare1 = function(design1.df, contr){
  n = nrow(design1.df)
  Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))
   
  
  pairMat = contr
  
  
  nBlk = nlevels(design1.df$Run)
  nPlot = nlevels(design1.df$Tag)
  nCag = nlevels(design1.df$Cag)
  
  nAni = nlevels(design1.df$Ani)
  Run.mat = with(design1.df, as.matrix(table(1:n, Run)))
  Tag.mat = with(design1.df, as.matrix(table(1:n, Tag)))
  Cag.mat = with(design1.df, as.matrix(table(1:n, Cag)))
  Ani.mat = with(design1.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
  Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))
  
  C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
  C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
  cage.Rep = n/nCag
  
  C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
  
  Pb = projMat(Run.mat)
  Pb1 = projMat(Tag.mat)
  
  blk.proj = (identityMat(n) - Pb) %*% (identityMat(n) - Pb1)
  
  info.mat = matMulti(blk.proj, Ani.mat, C.cage)
  
    if (any(abs(info.mat) > 1e-07)) {
        PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)),
            C.cage, Ani.mat)
  
        PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
    } else {
        PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani,
            Ani.mat)
  
    }
  
  contr.e = (apply(pairMat, 2, function(c)  t(c) %*% ginv(t(Trt.mat) %*% Trt.mat) %*% c)/apply(pairMat, 2, function(c)  t(c) %*% ginv(t(Trt.mat) %*% PP1 %*% Trt.mat) %*% c))
  
  
  return(contr.e)
}


#########################################################################################


round(basicContr(design16)[,1:7], 4)

contrCompare(design16, basicContr(design16)[,1:7])

1/mean(1/contrCompare(design16, basicContr(design16)[,1:7]))


round(basicContr(design14), 4)

contrCompare(design14, basicContr(design14)[,1:7])


1/mean(1/contrCompare(design14, basicContr(design14)[,1:7]))

ans14 = optim(ind, obj, swap, method = "SANN", design= design14, contr = basicContr(design14)[,1:7],
             control = list(maxit = 10000, temp = 0.1, trace = TRUE,
                            REPORT = 100))

contrCompare(design14,  asicContr(design14)[ans14$par,1:7])

ans16 = optim(ind, obj, swap, method = "SANN", design= design16, contr = basicContr(design16)[,1:7],
             control = list(maxit = 10000, temp = .1, trace = TRUE,
                            REPORT = 100))

contrCompare(design16,  basicContr(design16)[ans16$par,1:7])

######################################################################################## 
contr1 = c(1, -1, 0,0,0,0,0,0)
contr2 = c( 0,0,1, -1,0,0,0,0)
contr3 = c( 0,0,0,0,1, -1,0,0)
contr4 = c(0,0, 0,0,0,0,1, -1)
contr5 = c(1,1, -1,-1,0,0,0, 0)
contr6 = c(0,0, 0,0,1,1, -1,-1)
contr7 = c(1,1, 1, 1, -1, -1, -1, -1)
  
contr = cbind(contr1, contr2, contr3, contr4, contr5, contr6, contr7)/2  


contrCompare(design14, contr)

1/mean(1/contrCompare(design14, contr))

contrCompare(design16, contr)

1/mean(1/contrCompare(design16, contr))

obj = function(ind, design, contr){
 temp = contrCompare(design,  contr[ind,])
 
 sd(temp)
 #return(-max(temp))
 #return(-max(table(temp)))
}

swap = function(ind, design, contr) sample(length(ind),length(ind))

ans14 = optim(ind, obj, swap, method = "SANN", design= design14, contr = contr,
             control = list(maxit = 10000, temp = 0.1, trace = TRUE,
                            REPORT = 100))

contrCompare(design14,  contr[ans14$par,])

ans16 = optim(ind, obj, swap, method = "SANN", design= design16, contr = contr,
             control = list(maxit = 10000, temp = .1, trace = TRUE,
                            REPORT = 100))

contrCompare(design16,  contr[ans16$par,])



################################################################################
> contr = cbind(contr1, contr2, contr3, contr4, contr5, contr6, contr7)  
> contr
     contr1 contr2 contr3 contr4 contr5 contr6 contr7
[1,]      1      0      0      0      1      0      1
[2,]     -1      0      0      0      1      0      1
[3,]      0      1      0      0     -1      0      1
[4,]      0     -1      0      0     -1      0      1
[5,]      0      0      1      0      0      1     -1
[6,]      0      0     -1      0      0      1     -1
[7,]      0      0      0      1      0     -1     -1
[8,]      0      0      0     -1      0     -1     -1
> contrCompare(design14, contr)
   contr1    contr2    contr3    contr4    contr5    contr6    contr7 
0.8705128 0.8705128 0.8705128 0.8705128 0.7780749 0.7780749 0.9326923 
> fractions(contrCompare(design14, contr))
 contr1  contr2  contr3  contr4  contr5  contr6  contr7 
679/780 679/780 679/780 679/780 291/374 291/374  97/104 
> fractions(contrCompare(design16, contr))
contr1 contr2 contr3 contr4 contr5 contr6 contr7 
   6/7    4/5  12/13    6/7    8/9  24/31    4/5 
> 
> 1/mean(1/contrCompare(design14, contr))
[1] 0.8497616
> 1/mean(1/contrCompare(design16, contr))
[1] 0.84

################################################################################

contr1 =  c(1,1, 1, 1, -1, -1, -1, -1)
contr2 =  c(1,1, -1, -1, 1, 1, -1, -1)
contr3 =  c(1,-1, 1, -1, 1, -1,1, -1)
contr4 = contr1 * contr2
contr5 = contr1 * contr3
contr6 = contr2 * contr3
contr7 = contr1 * contr2 * contr3 
  
contr = cbind(contr1, contr2, contr3, contr4, contr5, contr6, contr7)  

contr
contrCompare(design14, contr)
contrCompare(design16, contr)

1/mean(1/contrCompare(design14, contr))
1/mean(1/contrCompare(design16, contr))


ans14 = optim(sample(ind), obj, swap, method = "SANN", design= design14, contr = contr,
             control = list(maxit = 10000, temp = 1, trace = TRUE,
                            REPORT = 100))

contrCompare(design14,  contr[ans14$par,])

ans16 = optim(sample(ind), obj, swap, method = "SANN", design= design16, contr = contr,
             control = list(maxit = 10000, temp = 1, trace = TRUE,
                            REPORT = 100))

contrCompare(design16,  contr[ans16$par,])






states = lapply(1:8, function(x) t(cbind(x, (x+1):8)))
states = states[-length(states)]

states= t(matrix(c(states, recursive=TRUE), nrow = 2))

contr.e = contrCompare(design16, pairwiseContr(design16))

contrastMat = matrix(0, ncol= 8, nrow = 8)

for( i in 1:length(contr.e)){   
  contrastMat[states[i,2], states[i,1]] = contr.e[i]
}

contrastMat


#Pick the best 4 pairwise contrasts 
contr1 = c(1, 0, 0,-1,0,0,0,0)
contr2 = c( 0,1,0, 0,0,0,-1,0)
contr3 = c( 0,0,1,0,0, 0,0,-1)
contr4 = c(0,0, 0,0,1,-1,0, 0)
contr5 = c(1,-1, 0, 1,0,0,-1, 0)
contr6 = c(0,0, 1,0,-1,-1, 0,1)
contr7 = c(1,1, -1, 1, -1,-1, 1, -1)
  
contr = cbind(contr1, contr2, contr3, contr4, contr5, contr6, contr7)/2  

contrCompare(design16, contr)

1/mean(1/contrCompare(design16, contr))


contr.e = contrCompare(design14, pairwiseContr(design14))

contrastMat = matrix(0, ncol= 8, nrow = 8)

for( i in 1:length(contr.e)){   
  contrastMat[states[i,2], states[i,1]] = contr.e[i]
}

contrastMat

contr1 = c(1, -1, 0,0,0,0,0,0)
contr2 = c( 0,0,1, -1,0,0,0,0)
contr3 = c( 0,0,0,0,1, -1,0,0)
contr4 = c(0,0, 0,0,0,0,1, -1)
contr5 = c(1,1, -1,-1,0,0,0, 0)
contr6 = c(0,0, 0,0,1,1, -1,-1)
contr7 = c(1,1, 1, 1, -1, -1, -1, -1)
  
contr = cbind(contr1, contr2, contr3, contr4, contr5, contr6, contr7)/2  


contrCompare(design14, contr)


1/mean(1/contrCompare(design14, contr))
