
library(infoDecompuTE)
# Load in the designs #####

load("designTag2.Rdata")

worstResDes <-  designTag1.list[[48789]]$design

midDes <-   designTag1.list[[207259]]$design

displayDes(midDes)

displayDes(worstResDes)

bestDes <- designTag1.list[[854022]]$design

bestDes$Plant[bestDes$Tray == 2] <-
  LETTERS[as.numeric(bestDes$Plant)[bestDes$Tray == 2] - 6]

summaryAovTwoPhase(
  bestDes,
  blk.str1 = "Tray/Plant",
  blk.str2 = "Run",
  trt.str = "Tag + Trt"
)$ANOVA



midDes$Plant[midDes$Tray == 2] <-
  LETTERS[as.numeric(midDes$Plant)[midDes$Tray == 2] - 6]

summaryAovTwoPhase(
  midDes,
  blk.str1 = "Tray/Plant",
  blk.str2 = "Run",
  trt.str = "Tag + Trt"
)$ANOVA


worstResDes$Plant[worstResDes$Tray == 2] <-
  LETTERS[as.numeric(worstResDes$Plant)[worstResDes$Tray == 2] - 6]


summaryAovTwoPhase(
  worstResDes,
  blk.str1 = "Tray/Plant",
  blk.str2 = "Run",
  trt.str = "Tag + Trt"
)$ANOVA


getVarTrtBest <- function(des, VCs, trts){
  
  des$Run <- factor(des$Run)
  des$Tag <- factor(des$Tag)
  des$Tray <- factor(des$Tray)
  des$Plant <- factor(des$Plant)
  des$Trt <- factor(des$Trt)
  
  n <- nrow(des)
  sigmaS <- VCs[4]
  sigmaP <- VCs[3]
  sigmaT <- VCs[2]
  sigmaR <- VCs[1]
  
  xTag <- c(1, 1, 1, 1)
  xTrt <- trts

  y <- sapply(1:1000, function(x) {
    res <- 10 + rnorm(nlevels(des$Run), sd = sigmaR)[des$Run] +
      rnorm(nlevels(des$Tray), sd = sigmaT)[des$Tray] +
      rnorm(nlevels(des$Plant), sd = sigmaP)[des$Plant] +
      rnorm(n, sd = sigmaS) +
      xTag[des$Tag] + xTrt[des$Trt]; 
    
    aov.Run <- aov(res ~ Tray / Plant + Error(Run), des);
    
    aov.Tray <- aov(res ~ Tag + Trt + Error(Run +  Tray / Plant), des);
    
    as.numeric(c(summary(aov.Run$Run)[[1]][2,3], 
      summary(aov.Tray$Tray)[[1]][3], 
      summary(aov.Tray$`Tray:Plant`)[[1]][2,3], 
      summary(aov.Tray$Within)[[1]][2,3], 
      summary.lm(aov.Tray$`Tray:Plant`)$coef[,1], 
      summary.lm(aov.Tray$`Tray:Plant`)$coef[,2]))
    })

    t(y)
}

getVarTrtMid <-function(des, VCs, trts){
  
  des$Run <- factor(des$Run)
  des$Tag <- factor(des$Tag)
  des$Tray <- factor(des$Tray)
  des$Plant <- factor(des$Plant)
  des$Trt <- factor(des$Trt)
  
  n <- nrow(des)
  sigmaS <- VCs[4]
  sigmaP <- VCs[3]
  sigmaT <- VCs[2]
  sigmaR <- VCs[1]
  
  xTag <- c(1, 1, 1, 1)
  xTrt <- trts

  y <- sapply(1:1000, function(x) {
    res <- 10 + rnorm(nlevels(des$Run), sd = sigmaR)[des$Run] +
      rnorm(nlevels(des$Tray), sd = sigmaT)[des$Tray] +
      rnorm(nlevels(des$Plant), sd = sigmaP)[des$Plant] +
      rnorm(n, sd = sigmaS) +
      xTag[des$Tag] + xTrt[des$Trt]; 
    
    aov.Run <- aov(res ~ Tray / Plant + Error(Run), des);
    
    aov.Tray <- aov(res ~ Tag + Trt + Error(Run +  Tray / Plant), des);
    
    as.numeric(c(summary(aov.Run$Run)[[1]][2,3], 
                 summary(aov.Tray$Tray)[[1]][3], 
                 summary(aov.Tray$`Tray:Plant`)[[1]][3,3], 
                 summary(aov.Tray$Within)[[1]][2,3], 
                 summary.lm(aov.Tray$`Tray:Plant`)$coef[,1], 
                 summary.lm(aov.Tray$`Tray:Plant`)$coef[,2]))
  })
  
  t(y)
}

getVarTrtWorst <- function(des, VCs, trts){
  
  des$Run <- factor(des$Run)
  des$Tag <- factor(des$Tag)
  des$Tray <- factor(des$Tray)
  des$Plant <- factor(des$Plant)
  des$Trt <- factor(des$Trt)
  
  n <- nrow(des)
  sigmaS <- VCs[4]
  sigmaP <- VCs[3]
  sigmaT <- VCs[2]
  sigmaR <- VCs[1]
  
  xTag <- c(1, 1, 1, 1)
  xTrt <- trts
  
  y <- sapply(1:1000, function(x) {
    res <- 10 + rnorm(nlevels(des$Run), sd = sigmaR)[des$Run] +
      rnorm(nlevels(des$Tray), sd = sigmaT)[des$Tray] +
      rnorm(nlevels(des$Plant), sd = sigmaP)[des$Plant] +
      rnorm(n, sd = sigmaS) +
      xTag[des$Tag] + xTrt[des$Trt]; 
    
    aov.Run <- aov(res ~ Tray / Plant + Error(Run), des);
    
    aov.Tray <- aov(res ~ Tag + Trt + Error(Run +  Tray / Plant), des);
    
    as.numeric(c(summary(aov.Run$Run)[[1]][2,3], 
                 summary(aov.Tray$Tray)[[1]][3], 
                 summary(aov.Tray$`Tray:Plant`)[[1]][3,3], 
                 summary(aov.Tray$Within)[[1]][2,3], 
                 summary.lm(aov.Tray$`Tray:Plant`)$coef[,1], 
                  summary.lm(aov.Tray$`Tray:Plant`)$coef[,2]))
  })
  
  t(y)
}

####################################################################################


#VCs <- c(0.0736, 0.0055, 0.0027, 0.0181)

VCs1 <- c(0.08, 0.006, 0.003, 0.02)
VCs2 <- c(0.08, 0.006, 0.003, 0.002)
VCs3 <- c(0.08, 0.06, 0.03, 0.02)
VCs4 <- c(0.08, 0.06, 0.003, 0.02)
VCs5 <- c(0.8, 0.006, 0.003, 0.02)

VCs <- VCs1

trts <- c(1,3,6)

compareVCs <- 
  function(VCs, trts){

  aovListBest <- getVarTrtBest(bestDes, VCs, trts)
  aovListMid <- getVarTrtMid(midDes, VCs, trts)
  aovListworst <- getVarTrtWorst(worstResDes, VCs, trts)
  
  ####################################################################
  
  
  runVCBest <- (aovListBest[,1] - aovListBest[,4])/4
  trayVCBest <- (aovListBest[,2] - aovListBest[,3])/12
  
  plantVCBest <- (aovListBest[,3] - aovListBest[,4])/2
  VCBest <- aovListBest[,4]
  
  
  runVCMid <- (aovListMid[,1] - aovListMid[,4])/4
  trayVCMid <- (aovListMid[,2] - aovListMid[,3])/12
  plantVCMid <- (aovListMid[,3] - aovListMid[,4])/2
  VCMid <- aovListMid[,4]
  
  runVCworst <- (aovListworst[,1] - aovListworst[,4])/4
  trayVCworst <- (aovListworst[,2] - aovListworst[,3])/12
  plantVCworst <- (aovListworst[,3] - aovListworst[,4])/2
  VCworst <- aovListworst[,4]
  
  
  Res <- 
  data.frame(
    Actual = VCs, 
    Best =
      c(sqrt(mean(runVCBest, na.rm = TRUE)),
        sqrt(mean(trayVCBest, na.rm = TRUE)),
        sqrt(mean(plantVCBest, na.rm = TRUE)),
        sqrt(mean(VCBest, na.rm = TRUE))),
    Mid =
      c(sqrt(mean(runVCMid, na.rm = TRUE)),
    sqrt(mean(trayVCMid, na.rm = TRUE)),
    sqrt(mean(plantVCMid, na.rm = TRUE)),
    sqrt(mean(VCMid, na.rm = TRUE))),
    Worst =
      c(sqrt(mean(runVCworst, na.rm = TRUE)),
    sqrt(mean(trayVCworst, na.rm = TRUE)),
    sqrt(mean(plantVCworst, na.rm = TRUE)),
    sqrt(mean(VCworst, na.rm = TRUE)))
  )
  
  rownames(Res) <- c("Run", "Tray", "Plant", "e")
  
  round(Res, 4)
}

set.seed(08032018)
overallRes <- 
  rbind(
    compareVCs(VCs1, trts), 
    compareVCs(VCs2, trts), 
    compareVCs(VCs3, trts), 
    compareVCs(VCs4, trts), 
    compareVCs(VCs5, trts))


xtable::xtable(overallRes, digits = 4)

####################################################################################

compareTrts <- 
  function(VCs, trts){
    
    aovListBest <- getVarTrtBest(bestDes, VCs, trts)
    aovListMid <- getVarTrtMid(midDes, VCs, trts)
    aovListworst <- getVarTrtWorst(worstResDes, VCs, trts)
    
    ####################################################################
  
    Res <- 
      data.frame(
        Actual = c(trts[2]-trts[1], trts[3]-trts[1]), 
        Best =
          c(mean(aovListBest[,5]), mean(aovListBest[,6])),
        BestSE = c(mean(aovListBest[,7]), mean(aovListBest[,8])),
        Mid =
          c(mean(aovListMid[,6]), mean(aovListMid[,7])),
        MidSE = c(mean(aovListMid[,9]), mean(aovListMid[,10])),
        Worst =
          c(mean(aovListworst[,7]), mean(aovListworst[,8])), 
        WorstSE = c(mean(aovListworst[,11]), mean(aovListworst[,12]))
      )
    
    rownames(Res) <- c("Trt2 - Trt1", "Trt3 - Trt1")
    
    round(Res, 4)
  }

trts1 <- c(1,3,6)
trts2 <- c(1,2,4)
trts3 <- c(6,3,1)
trts4 <- c(4,2,1)
trts5 <- c(1,1,1)

set.seed(08032018)
overallResTrt <- 
  rbind(
    compareTrts(VCs1, trts1), 
    compareTrts(VCs1, trts2), 
    compareTrts(VCs1, trts3), 
    compareTrts(VCs1, trts4), 
    compareTrts(VCs1, trts5))

xtable::xtable(overallResTrt, digits = 4)


boxplot(data.frame(Best = aovListBest[,3], 
                   Mid = aovListMid[,3], 
                   Worst = aovListworst[,3]), 
        main = "Residual Mean square")


trts <- c(1,3,6)

set.seed(08032018)
aovListBest <- getVarTrtBest(bestDes, VCs, trts)
aovListMid <- getVarTrtMid(midDes, VCs, trts)
aovListworst <- getVarTrtWorst(worstResDes, VCs, trts)

runVCBest <- (aovListBest[,1] - aovListBest[,4])/4
trayVCBest <- (aovListBest[,2] - aovListBest[,3])/12

plantVCBest <- (aovListBest[,3] - aovListBest[,4])/2
VCBest <- aovListBest[,4]


runVCMid <- (aovListMid[,1] - aovListMid[,4])/4
trayVCMid <- (aovListMid[,2] - aovListMid[,3])/12
plantVCMid <- (aovListMid[,3] - aovListMid[,4])/2
VCMid <- aovListMid[,4]

runVCworst <- (aovListworst[,1] - aovListworst[,4])/4
trayVCworst <- (aovListworst[,2] - aovListworst[,3])/12
plantVCworst <- (aovListworst[,3] - aovListworst[,4])/2
VCworst <- aovListworst[,4]


boxplot(data.frame(Best = aovListBest[,3], 
                   Mid = aovListMid[,3], 
                   Worst = aovListworst[,3]), 
        main = "Residual Mean square")



# SED ####

opar <-par(mfrow = c(2,2))

boxplot(data.frame(Best = sqrt(runVCBest), 
                   Mid = sqrt(runVCMid), 
                   Worst = sqrt(runVCworst)), 
        main = "VC estimates for Runs")
#abline(h = VCs[1], col = "red")


boxplot(data.frame(Best = sqrt(trayVCBest), 
                   Mid = sqrt(trayVCMid), 
                   Worst = sqrt(trayVCworst)), 
        main = "VC estimates for Trays")
#abline(h = VCs[2], col = "red")

boxplot(data.frame(Best = sqrt(plantVCBest), 
                   Mid = sqrt(plantVCMid), 
                   Worst = sqrt(plantVCworst)), 
        main = "VC estimates for Plants")
#abline(h = VCs[3], col = "red")


boxplot(data.frame(Best = sqrt(VCBest), 
                   Mid = sqrt(VCMid), 
                   Worst = sqrt(VCworst)), 
        main = "VC estimates for Measurment Errors")
#abline(h = VCs[4], col = "red")
par(opar)




sd(sqrt(runVCBest), na.rm = TRUE)
sd(sqrt(runVCMid), na.rm = TRUE)
sd(sqrt(runVCworst), na.rm = TRUE)


sd(sqrt(trayVCBest), na.rm = TRUE)
sd(sqrt(trayVCMid), na.rm = TRUE)
sd(sqrt(trayVCworst), na.rm = TRUE)

sd(sqrt(plantVCBest), na.rm = TRUE)
sd(sqrt(plantVCMid), na.rm = TRUE)
sd(sqrt(plantVCworst), na.rm = TRUE)

sd(sqrt(VCBest), na.rm = TRUE)
sd(sqrt(VCMid), na.rm = TRUE)
sd(sqrt(VCworst), na.rm = TRUE)




boxplot(data.frame(Best = aovListBest[,5], 
                   Mid = aovListMid[,6], 
                   Worst = aovListworst[,7]), 
        main = "VC estimates for Measurment Errors")

boxplot(data.frame(Best = aovListBest[,6], 
                   Mid = aovListMid[,7], 
                   Worst = aovListworst[,8]), 
        main = "VC estimates for Measurment Errors")














