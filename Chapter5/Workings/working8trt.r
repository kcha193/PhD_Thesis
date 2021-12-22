

 sourceDir <- function(path, trace = TRUE, ...) {
 library(MASS)
library(lattice)
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
       if(trace) cat(nm,":")           
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
 }


sourceDir(path = "C:/Users/kcha193/My Dropbox/R functions for two-phase experiments")



summary.aov.twoPhase(design.df = phase2designEX5, blk.str1 = "Ani", blk.str2 = "Run",
trt.str = "Trt + Tag")

summary.aov.modified(response = y, blk.str2 = "Set/Array", blk.str1 = "Plant",
                 trt.str = "Treat + Dye", design = design.df)


asreml.fit = asreml(y ~ Trt + Tag, random = ~ Run + Ani, data = design)

summary(asreml.fit)
svc.asreml(asreml.fit)



design.df = data.frame(Set = factor(rep(1:4, each = 4)),
                       Array = factor(rep(1:2, time = 4, each = 2)),
                       Dye = factor(c("Gre", "Red")),
                       Treat = factor(LETTERS[c(1,2,2,1,
                                                1,2,2,1,
                                                1,2,2,1,
                                                1,2,2,1)]),
                       Plant = interaction(factor(rep(1:4, each = 4)),
                                           factor(LETTERS[c(1,2,2,1,
                                                            1,2,2,1,
                                                            1,2,2,1,
                                                            1,2,2,1)])))


aov.table = summary.aov.twoPhase(design.df , blk.str2 = "Set/Array", blk.str1 = "Plant", 
                trt.str = "Treat + Dye", var.comp = c("Plant", "Set:Array"),response = y)

getVcEDF(aov.table = aov.table, row.MS = row.MS, true.VC = true.VC, neg.VC = neg.VC)

summary.aov.onePhase(design.df , blk.str = "Set/Array", trt.str  = "Plant", 
                var.comp = c("Set:Array"), response = rnorm(16))

sim9 = simData1(design =  design.df, 
            blk.str2 = "Set/Array", blk.str1 = "Plant", 
            trt.str = "Treat + Dye", var.comp = c("Plant", "Set:Array"), 
            gamma.A = c(0, 0.25, 1, 4, 100),          
            gamma.P =c(10^((-8:8)/2)), 
            VC.resid = 0.1,
            nS = 4, nVc = 3, nSim = 10000, neg.VC = FALSE)                        

sim9 = simData1(design =  design.df, 
            blk.str2 = "Set/Array", blk.str1 = "Plant", 
            trt.str = "Treat + Dye", var.comp = c("Plant", "Set:Array"), 
            gamma.A = c(0, 0.25),          
            gamma.P =c(10^((-8:8)/2)), 
            VC.resid = 0.1,
            nS = 4, nVc = 3, nSim = 1000, neg.VC = FALSE)                        

index = with(sim9$t, paste("gamma.array = ", gamma.array, ", gamma.plant = ", gamma.plant, sep = ""))


colmean = apply(sim9$t[,-c(1:2)], 2, function(x) tapply(x, index, function(y) round(mean(y), 5)))
apply(sim9$t[,-c(1:2)], 2, function(x) tapply(x, index, length))

plotEDF(sim9)

plotEDF1(sim5)

backup.sim5 = sim5

REML.VC = matrix(sim5$tempV$REML.VC, nrow = 3)
REML.VC[which(REML.VC<0)]=0
LC.VC = matrix(sim5$tempV$LC.VC, nrow = 3)
LC.VC[which(LC.VC<0)]=0
REAL.VC = matrix(sim5$tempV$TRUE.VC, nrow = 3)

sim5$tempS$REML.EDF = as.numeric(apply(REML.VC, 2, function(x) 
  getVcEDF(aov.table = aov.table, row.MS = NA, true.VC = x, neg.VC = TRUE)$Stratum[,"TRUE.EDF"]))
sim5$tempS$LC.EDF = as.numeric(apply(LC.VC, 2, function(x) 
  getVcEDF(aov.table = aov.table, row.MS = NA, true.VC = x, neg.VC = TRUE)$Stratum[,"TRUE.EDF"]))

plotEDF1(sim5)

################################################################################


 names(sim5)
 [1] "gamma.array" "gamma.plant" "e"           "Plant"       "Set:Array"   "e"          
 [7] "Plant"       "Set:Array"   "e"           "Plant"       "Set:Array"  


index = with(sim5, paste("gamma.array = ", gamma.array, ", gamma.plant = ", gamma.plant, sep = ""))

unique(index)
histVCplots(sim5, index, 1)

histVCplots(sim3, index, 1)
histVCplots(sim3, index, 10)
histVCplots(sim3, index, 19)


######################################################################################################
gamma.array =  1 
gamma.plant =  1 
> apply(sim7[,-c(1:2)], 2, length)
          e       Plant   Set:Array         e.1     Plant.1 Set:Array.1         e.2 
       9976        9976        9976        9976        9976        9976        9976 
    Plant.2 Set:Array.2 
       9976        9976 
> apply(sim7[,-c(1:2)], 2, mean)
          e       Plant   Set:Array         e.1     Plant.1 Set:Array.1         e.2 
  1.0969961   0.9999697   0.9650335   0.9871548   1.0325913   1.0094809   1.0000000 
    Plant.2 Set:Array.2 
  1.0000000   1.0000000 



> apply(sim7[,-c(1:2)], 2, length)
          e       Plant   Set:Array         e.1     Plant.1 Set:Array.1         e.2 
      99762       99762       99762       99762       99762       99762       99762 
    Plant.2 Set:Array.2 
      99762       99762 
> apply(sim7[,-c(1:2)], 2, mean)
          e       Plant   Set:Array         e.1     Plant.1 Set:Array.1         e.2 
  1.1107537   0.9806242   0.9665426   0.9981220   1.0055178   1.0031882   1.0000000 
    Plant.2 Set:Array.2 
  1.0000000   1.0000000 

gamma.array =  1 
gamma.plant =  10 
> apply(sim7[,-c(1:2)], 2, length)
          e       Plant   Set:Array         e.1     Plant.1 Set:Array.1         e.2 
       9999        9999        9999        9999        9999        9999        9999 
    Plant.2 Set:Array.2 
       9999        9999 
> apply(sim7[,-c(1:2)], 2, mean)
          e       Plant   Set:Array         e.1     Plant.1 Set:Array.1         e.2 
  1.0279407   9.9827132   1.0484055   1.0070925   9.9474858   0.9819257   1.0000000 
    Plant.2 Set:Array.2 
 10.0000000   1.0000000
 
 
 
gamma.array =  0.25 
gamma.plant =  0.001 
VC.resid = 1

> apply(sim7$t[,-c(1:2)], 2, length)
          e       Plant   Set:Array         e.1     Plant.1 Set:Array.1         e.2 
       9841        9841        9841        9841        9841        9841        9841 
    Plant.2 Set:Array.2 
       9841        9841 
> apply(sim7$t[,-c(1:2)], 2, mean)
          e       Plant   Set:Array         e.1     Plant.1 Set:Array.1         e.2 
1.050153221 0.003211987 0.238151550 0.999037687 0.002324175 0.263253312 1.000000000 
    Plant.2 Set:Array.2 
0.001000000 0.250000000 

gamma.array =  0.25 
gamma.plant =  0.001 
VC.resid = 0.001

> apply(sim7$t[,-c(1:2)], 2, length)
          e       Plant   Set:Array         e.1     Plant.1 Set:Array.1 
       9845        9845        9845        9845        9845        9845 
        e.2     Plant.2 Set:Array.2 
       9845        9845        9845 
> apply(sim7$t[,-c(1:2)], 2, mean)
           e        Plant    Set:Array          e.1      Plant.1  Set:Array.1 
1.045293e-03 4.772950e-06 2.386972e-04 9.996445e-04 6.518537e-06 2.613134e-04 
         e.2      Plant.2  Set:Array.2 
1.000000e-03 1.000000e-06 2.500000e-04 



 nVc = 3
REML1.VC = numeric(nVc); REML2.VC = numeric(nVc); LC.VC   = numeric(nVc); TRUE.VC  = numeric(nVc);

#Fail t converge
for( i in 1: (nrow(sim10$y)-1)){
print(i)
  aov.table = summary.aov.twoPhase(design.df , blk.str2 = "Set/Array", blk.str1 = "Plant", 
                  trt.str = "Treat + Dye", var.comp = c("Plant", "Set:Array"),response = sim10$y[i,])
  
  tmp = getFailVc(aov.table = aov.table, row.MS = NA, 
  true.VC =  c(VC.resid, (gamma.plant * VC.resid),(gamma.array * VC.resid)),
   neg.VC =FALSE)
   
   tmpV = tmp$V
   REML1.VC   =  rbind(REML1.VC,      tmpV[,"REML.var.comp1"])
   REML2.VC   =  rbind(REML2.VC,      tmpV[,"REML.var.comp2"])
   LC.VC     =  rbind(LC.VC  ,      tmpV[,"LC.var.comp"]    )
   TRUE.VC   =  rbind(TRUE.VC,      tmpV[,"TRUE.var.comp"]  ) 
}

REML1.VC   =  REML1.VC[-1,]; REML2.VC   =  REML2.VC[-1,]; LC.VC   =  LC.VC[-1,]; TRUE.VC   = TRUE.VC[-1,] 

apply(REML1.VC, 2, mean)
apply(REML2.VC, 2, mean)
apply(rbind( REML1.VC, REML2.VC), 2, mean)

apply(LC.VC, 2, mean)

apply(rbind(sim7$t[,c(3,4,5)], REML1.VC, REML2.VC), 2, mean)
apply(rbind(sim7$t[,c(3,4,5)], REML1.VC), 2, mean)
apply(rbind(sim7$t[,c(3,4,5)], REML2.VC), 2, mean)

apply(rbind(sim7$t[,c(3,4,5)],(REML1.VC + REML2.VC)/2), 2, mean)

apply(rbind(sim7$t[,c(6,7,8)], LC.VC), 2, mean)

gamma.array =  0 
gamma.plant =  0.001
VC.resid = 0.001

> apply(sim7$t[,-c(1:2)], 2, length)
          e       Plant   Set:Array         e.1     Plant.1 Set:Array.1 
       9846        9846        9846        9846        9846        9846 
        e.2     Plant.2 Set:Array.2 
       9846        9846        9846 
> apply(sim7$t[,-c(1:2)], 2, mean)
            e         Plant     Set:Array           e.1       Plant.1 
 1.023993e-03  1.086048e-05 -6.470956e-06  9.922405e-04  1.013102e-05 
  Set:Array.1           e.2       Plant.2   Set:Array.2 
 3.624314e-06  1.000000e-03  1.000000e-06  0.000000e+00 

gamma.array =  0 
gamma.plant =  0.001
VC.resid = 100
> apply(sim7$t[,-c(1:2)], 2, length)
          e       Plant   Set:Array         e.1     Plant.1 Set:Array.1 
       9848        9848        9848        9848        9848        9848 
        e.2     Plant.2 Set:Array.2 
       9848        9848        9848 
> apply(sim7$t[,-c(1:2)], 2, mean)
          e       Plant   Set:Array         e.1     Plant.1 Set:Array.1 
104.1305167   0.1350152  -1.3403010 100.8560742   0.2796711  -0.7643562 
        e.2     Plant.2 Set:Array.2 
100.0000000   0.1000000   0.0000000 

################################################################################
   
sim2 = simData(design =  phase2designEX5, 
            blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Trt + Tag", 
            gamma.R = c(0, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000), 
            gamma.A =c(10^((-10:10)/2)), 
            VC.resid = 0.001,
            nS = 4, nVc = 3, nSim = 10000, neg.VC = TRUE)                        

plotEDF(sim2)







