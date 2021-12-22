

library(ggplot2)
library(xtable)

load("Working.Rdata")

names(CRD232Tag4.EDF$tempS)[3:4] = c("REML", "LC")

edfPlot(CRD232Tag4.EDF)  + theme(legend.position = "top") + ggtitle("")
ggsave("CRD232.pdf", width=12, height = 7) 



CRD262Tag4.des

summaryAovTwoPhase(
  design.df = CRD262Tag4.des[-c(1),],
  blk.str2 = "Run",
  blk.str1 = "Ani",
  trt.str = "Tag + Trt",
  decimal = TRUE,
  digits = 4
)


CRD262Tag4.missing = simDataChisq(design =  CRD262Tag4.des,
                          blk.str2 = "Run", blk.str1 = "Ani",
                          trt.str = "Tag + Trt",
                          gamma.R = c(0, 0.25, 1, 4, 100),
                          gamma.A = c(10^((-8:8)/2)),
                          VC.resid = 1,
                          nS = 4, nVc = 3, nSim = 100, seed = 1, neg.VC = TRUE)

edfPlot(CRD262Tag4.missing)  + theme(legend.position = "top") + ggtitle("")



edfPlot(CRD232Tag4.EDF)  + theme(legend.position = "top") + ggtitle("")
ggsave("CRD232.pdf", width=12, height = 7) 



######

CRD262Tag4.des

summaryAovTwoPhase(
  design.df = CRD262Tag4.des,
  blk.str2 = "Run",
  blk.str1 = "Ani",
  trt.str = "Tag + Trt",
  decimal = TRUE,
  digits = 4
)

design.summary.CRD(CRD262Tag8.des)

  n = nrow(CRD262Tag8.des)
Run.mat = with(CRD262Tag8.des, as.matrix(table(1:n, Run)))
Tag.mat = with(CRD262Tag8.des, as.matrix(table(1:n, Tag)))
Ani.mat = with(CRD262Tag8.des, as.matrix(table(1:n, Ani)))
Trt.mat = with(CRD262Tag8.des, as.matrix(table(1:n, Trt)))

 with(design.df, table(Trt, Tag))

 t(Tag.mat) %*% Trt.mat
 

out <- matrix(paste0(CRD262Tag8.des$Ani, CRD262Tag8.des$Trt), ncol = nlevels(factor(CRD262Tag8.des$Tag)), 
              nrow = nlevels(CRD262Tag8.des$Run), byrow = TRUE)
out

xtable(out)

summaryAovTwoPhase(
  design.df = CRD262Tag8.des,
  blk.str2 = "Run",
  blk.str1 = "Ani",
  trt.str = "Tag + Trt",
  decimal = TRUE,
  digits = 4
)

out <- matrix(paste0(CRD262Tag4.des$Ani, CRD262Tag4.des$Trt), ncol = nlevels(factor(CRD262Tag4.des$Tag)), 
              nrow = nlevels(CRD262Tag4.des$Run), byrow = TRUE)
out

xtable(out)

###### 

colnames(CRD262Tag4$tempS) <- gsub(".EDF", "", colnames(CRD262Tag4$tempS) )
colnames(CRD262Tag8$tempS) <- gsub(".EDF", "", colnames(CRD262Tag8$tempS) )

edfPlotCompare.REML(CRD262Tag4, CRD262Tag8,  compare = c("4-plex", "8-plex")) + 
  theme(legend.position = "top") + ggtitle("")
ggsave("CRD262.pdf", width=12, height = 7)  


#####

out <- matrix(paste0(CRD822Tag4.des$Ani, CRD822Tag4.des$Trt), ncol = nlevels(factor(CRD822Tag4.des$Tag)), 
              nrow = nlevels(CRD822Tag4.des$Run), byrow = TRUE)
out

xtable(out)

out <- matrix(paste0(CRD822Tag8.des$Ani, CRD822Tag8.des$Trt), ncol = nlevels(factor(CRD822Tag8.des$Tag)), 
              nrow = nlevels(CRD822Tag8.des$Run), byrow = TRUE)
out

xtable(out)

colnames(CRD822Tag4$tempS) <- gsub(".EDF", "", colnames(CRD822Tag4$tempS) )
colnames(CRD822Tag8$tempS) <- gsub(".EDF", "", colnames(CRD822Tag8$tempS) )

edfPlotCompare.REML(CRD822Tag4, CRD822Tag8,  compare = c("4-plex", "8-plex")) + 
  theme(legend.position = "top") + ggtitle("")
ggsave("CRD822.pdf", width=12, height = 7)  

############

out <- matrix(paste0(CRD462Tag4.des$Ani, CRD462Tag4.des$Trt), ncol = nlevels(factor(CRD462Tag4.des$Tag)), 
              nrow = nlevels(CRD462Tag4.des$Run), byrow = TRUE)
out

xtable(out)

out <- matrix(paste0(CRD462Tag8.des$Ani, CRD462Tag8.des$Trt), ncol = nlevels(factor(CRD462Tag8.des$Tag)), 
              nrow = nlevels(CRD462Tag8.des$Run), byrow = TRUE)
out

xtable(out)

colnames(CRD462Tag4$tempS) <- gsub(".EDF", "", colnames(CRD462Tag4$tempS) )
colnames(CRD462Tag8$tempS) <- gsub(".EDF", "", colnames(CRD462Tag8$tempS) )

edfPlotCompare.REML(CRD462Tag4, CRD462Tag8,  compare = c("4-plex", "8-plex")) + 
  theme(legend.position = "top") + ggtitle("")
ggsave("CRD462.pdf", width=12, height = 7)  


###############################################################################################
#Block stuff
#Compare 4 and 8 tags experiments RBD

RBD44424Run= optRBD(nTrt = 4, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 4, iter = 1000,
                   confoundCag = FALSE, upperValue = 1)


RBD44424Tag = optRBD(nTrt = 4, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 4, iter = 1000,
                    confoundCag = TRUE, upperValue = 1)

design.summary.RBD(RBD44424Run, FALSE)
design.summary.RBD(RBD44424Tag, FALSE)

RBD44428Run = optRBD(nTrt = 4, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 8, iter = 100,
                    confoundCag = FALSE, upperValue = 1)

RBD44428Tag = optRBD(nTrt = 4, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 8, iter = 1000,
                    confoundCag = TRUE, upperValue = 1)

design.summary.RBD(RBD44428Run, FALSE)
design.summary.RBD(RBD44428Tag, FALSE)

###################################################################################################

(RBD44424Run$ANI <- LETTERS[interaction(RBD44424Run$Ani, RBD44424Run$Cag)])

sort(unique(paste0(RBD44424Run$Cag, RBD44424Run$ANI)))

out <- matrix(paste0(RBD44424Run$Cag, RBD44424Run$ANI, RBD44424Run$Trt), ncol = nlevels(factor(RBD44424Run$Tag)), 
              nrow = nlevels(RBD44424Run$Run), byrow = TRUE)

out

xtable(out)

design.summary.RBD(RBD44424Run, FALSE)

RBD44424Run.EDF = simDataChisqPRBD(design =  RBD44424Run,
                             blk.str2 = "Run", blk.str1 = "Cag/Ani",
                             trt.str = "Tag + Trt",
                             gamma.R = c(0, 0.25, 1, 4, 100),
                             gamma.C = c(0,  0.25, 1, 4, 100),
                             gamma.A =c(10^((-8:8)/2)),
                             VC.resid = 1,
                             nS = 4, nVc = 4, nSim = 1000, seed = 1, neg.VC = TRUE, checkFun = FALSE)

(aov.table =summaryAovTwoPhase(RBD44424Run,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

( rowNames = extractName(rownames(aov.table)))

RBD44424Run.NA1 = simDataChisqPRBD(design =  RBD44424Run,
                             blk.str2 = "Run", blk.str1 = "Cag/Ani",
                             trt.str = "Tag + Trt",
                             gamma.R = c(NA),
                             gamma.C = NA,
                             gamma.A =c(10^((-8:8)/2)),
                             VC.resid = 1,
                             nS = 2, nVc = 2, nSim = 100, seed = 1, neg.VC = TRUE, checkFun = FALSE,
                 row.MS =rowNames[c(8, 11)])

RBD44424Run.NA2 = simDataChisqPRBD(design =  RBD44424Run,
                             blk.str2 = "Run", blk.str1 = "Cag/Ani",
                             trt.str = "Tag + Trt",
                             gamma.R = c(0,  0.25, 1, 4, 100),
                             gamma.C = NA,
                             gamma.A =c(10^((-8:8)/2)),
                             VC.resid = 1,
                             nS = 3, nVc = 2, nSim = 10, seed = 1, neg.VC = TRUE, checkFun = FALSE,
                 row.MS =rowNames[c(3, 8, 11)])

test <- RBD44424Run.NA2$tempS
test$gamma.cag <- test$gamma.run
test$gamma.run <- NA

temp <- rbind(RBD44424Run.EDF$tempS, RBD44424Run.NA2$tempS, test, RBD44424Run.NA1$tempS)




(RBD44424Tag$ANI <- LETTERS[interaction(RBD44424Tag$Ani, RBD44424Tag$Cag)])

sort(unique(paste0(RBD44424Tag$Cag, RBD44424Tag$ANI)))

out <- matrix(paste0(RBD44424Tag$Cag, RBD44424Tag$ANI, RBD44424Tag$Trt), 
              ncol = nlevels(factor(RBD44424Tag$Tag)), 
              nrow = nlevels(RBD44424Tag$Run), byrow = TRUE)

out

xtable(out)

RBD44424Tag.EDF = simDataChisqPRBD(design =  RBD44424Tag,
           blk.str2 = "Run", blk.str1 = "Cag/Ani",
           trt.str = "Tag + Trt",
           gamma.R = c(0, 0.25, 1, 4, 100),
           gamma.C = c(0,  0.25, 1, 4, 100),
           gamma.A =c(10^((-8:8)/2)),
           VC.resid = 1,
           nS = 6, nVc = 4, nSim = 1000, seed = 1, neg.VC = TRUE, checkFun = FALSE)

(aov.table =summaryAovTwoPhase(RBD44424Tag,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

( rowNames = extractName(rownames(aov.table)))


RBD44424Tag.NA = simDataChisqPRBD(design =  RBD44424Tag,
           blk.str2 = "Run", blk.str1 = "Cag/Ani",
           trt.str = "Tag + Trt",
           gamma.R = c(0, 0.25, 1, 4, 100),
           gamma.C = c(NA),
           gamma.A =c(10^((-8:8)/2)),
           VC.resid = 1,
           nS = 4, nVc = 3, nSim = 100, seed = 1, neg.VC = TRUE, checkFun = FALSE,
            row.MS = rowNames[c(3, 4, 11, 14)])

RBD44424Tag.NA1 = simDataChisqPRBD(design =  RBD44424Tag,
           blk.str2 = "Run", blk.str1 = "Cag/Ani",
           trt.str = "Tag + Trt",
           gamma.R = c(NA),
           gamma.C = c(0, 0.25, 1, 4, 100), 
           gamma.A =c(10^((-8:8)/2)),
           VC.resid = 1,
           nS = 3, nVc = 3, nSim = 100, seed = 1, neg.VC = TRUE, checkFun = FALSE,
            row.MS =rowNames[c(8, 11, 14)])

RBD44424Tag.NA2 = simDataChisqPRBD(design =  RBD44424Tag,
           blk.str2 = "Run", blk.str1 = "Cag/Ani",
           trt.str = "Tag + Trt",
           gamma.R = c(NA),
           gamma.C = c(NA), 
           gamma.A =c(10^((-8:8)/2)),
           VC.resid = 1,
           nS = 2, nVc = 2, nSim = 100, seed = 1, neg.VC = TRUE, checkFun = FALSE,
            row.MS =rowNames[c(11, 14)])

temp1 <- rbind(RBD44424Tag.EDF$tempS, RBD44424Tag.NA$tempS,
               RBD44424Tag.NA1$tempS, RBD44424Tag.NA2$tempS)


colnames(temp) <- gsub(".EDF", "", colnames(temp) )
colnames(temp1) <- gsub(".EDF", "", colnames(temp1) )

edfPlotCompareCag.REML(temp, temp1,  
                       compare = c("Confounded With Run", 
                                   "Confounded With Tag"))+ 
  theme(legend.position = "top") + ggtitle("")
ggsave("CRD44424.pdf", width=20, height = 12)  

###################################################################################################



(RBD44428Run$ANI <- LETTERS[interaction(RBD44428Run$Ani, RBD44428Run$Cag)])

sort(unique(paste0(RBD44428Run$Cag, RBD44428Run$ANI)))

out <- matrix(paste0(RBD44428Run$Cag, RBD44428Run$ANI, RBD44428Run$Trt), 
              ncol = nlevels(factor(RBD44428Run$Tag)), 
              nrow = nlevels(RBD44428Run$Run), byrow = TRUE)

out

xtable(out)

RBD44428Run.EDF = simDataChisqPRBD(design =  RBD44428Run,
                                   blk.str2 = "Run", blk.str1 = "Cag/Ani",
                                   trt.str = "Tag + Trt",
                                   gamma.R = c(0, 0.25, 1, 4, 100),
                                   gamma.C = c(0,  0.25, 1, 4, 100),
                                   gamma.A =c(10^((-8:8)/2)),
                                   VC.resid = 1,
                                   nS = 5, nVc = 4, nSim = 1000, seed = 1,
                                   neg.VC = TRUE, checkFun = FALSE)


temp3 <- RBD44428Run.EDF$tempS
temp31 <- RBD44428Run.EDF$tempS
temp32 <- RBD44428Run.EDF$tempS
temp33 <- RBD44428Run.EDF$tempS

temp31$gamma.run <- temp33$gamma.run <-NA
temp32$gamma.cag <- temp33$gamma.cag <-NA

temp3 <- rbind(unique(temp33), unique(temp32), unique(temp31), temp3)


(RBD44428Tag$ANI <- LETTERS[interaction(RBD44428Tag$Ani, RBD44428Tag$Cag)])

sort(unique(paste0(RBD44428Tag$Cag, RBD44428Tag$ANI)))

out <- matrix(paste0(RBD44428Tag$Cag, RBD44428Tag$ANI, RBD44428Tag$Trt), ncol = nlevels(factor(RBD44428Tag$Tag)), 
              nrow = nlevels(RBD44428Tag$Run), byrow = TRUE)

out

xtable(out)

RBD44428Tag.EDF = simDataChisqPRBD(design =  RBD44428Tag,
                                   blk.str2 = "Run", blk.str1 = "Cag/Ani",
                                   trt.str = "Tag + Trt",
                                   gamma.R = c(0, 0.25, 1, 4, 100),
                                   gamma.C = c(NA),
                                   gamma.A =c(10^((-8:8)/2)),
                                   VC.resid = 1,
                                   nS = 4, nVc = 3, nSim = 1000, seed = 1, 
                                   neg.VC = TRUE, checkFun = FALSE)


(aov.table =summaryAovTwoPhase(RBD44424Tag,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

( rowNames = extractName(rownames(aov.table)))


RBD44428TagNA.EDF = simDataChisqPRBD(design =  RBD44428Tag,
                                   blk.str2 = "Run", blk.str1 = "Cag/Ani",
                                   trt.str = "Tag + Trt",
                                   gamma.R = c(NA),
                                   gamma.C = c(NA),
                                   gamma.A =c(10^((-8:8)/2)),
                                   VC.resid = 1,
                                   nS = 2, nVc = 2, nSim = 10, seed = 1, 
                                   neg.VC = TRUE, checkFun = FALSE,
            row.MS =rowNames[c(11, 14)])

colnames(RBD44428TagNA.EDF$tempS) <- gsub(".EDF", "", colnames(RBD44428TagNA.EDF$tempS) )

temp4 <- RBD44428Run.EDF$tempS
temp4$REML <- NA
temp4$LC <- NA
temp41 <- temp31
temp41$REML <- NA
temp41$LC <- NA

temp4 <- rbind(RBD44428TagNA.EDF$tempS, RBD44428Tag.EDF$tempS, unique(temp41), temp4)

colnames(temp3) <- gsub(".EDF", "", colnames(temp3) )
colnames(RBD44428Tag.EDF$tempS) <- gsub(".EDF", "", colnames(RBD44428Tag.EDF$tempS) )

edfPlotCompareCag.REML(temp3, temp4,  
    compare = c("Confounded With Run", "Confounded With Tag"))+ 
  theme(legend.position = "top") + ggtitle("")
ggsave("CRD44428.pdf", width=20, height = 12)  

save.image("working0510.Rdata")



#####################################################################



edfPlotCompare1(temp, temp1,  temp3, temp4,
                       compare = c("Confounded With Run/4-plex", "Confounded With Tag/4-plex", 
                                   "Confounded With Run/8-plex", "Confounded With Tag/8-plex"))+ 
  theme(legend.position = "top") + ggtitle("")
ggsave("CRD44424vs8.pdf", width=12, height = 7.5)  



save.image("working0510.Rdata")



#####################################################################

load("working0510.Rdata")

library(optimTE)


bibd <- optBIBD(6,6,5,2,4,TRUE)

bibd65624 <- optBIBD(nTrt = 6, bRep = 5, nCag = 6, tRep = 2, nPlot = 4,
                    iter = 3000, confoundCag = FALSE)


(aov.table =summaryAovTwoPhase(RBD44428Tag[-(1:8),],  blk.str2 = "Run", blk.str1 = "Cag/Ani",
                               trt.str = "Tag + Trt")$A)

( rowNames = extractName(rownames(aov.table)))


bibd65624.EDF = simDataChisqPRBD(design =  RBD44428Tag[-(1:8),],
                                   blk.str2 = "Run", blk.str1 = "Cag/Ani",
                                   trt.str = "Tag + Trt",
                                   gamma.R = c(0, 0.25, 1, 4, 100),
                                   gamma.C = c(0,  0.25, 1, 4, 100),
                                   gamma.A =c(10^((-8:8)/2)),
                                   VC.resid = 1,
                                   nS = 5, nVc = 4, nSim = 10, seed = 1,
                                   neg.VC = TRUE, checkFun = FALSE)



edfPlotCompareCag.REML(bibd65624.EDF$tempS, bibd65624.EDF$tempS,  
                       compare = c("Confounded With Run", 
                                   "Confounded With Tag"))+ 
  theme(legend.position = "top") + ggtitle("")




