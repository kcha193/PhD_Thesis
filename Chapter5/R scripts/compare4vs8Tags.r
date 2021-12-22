#22/05/2013 02:26:39 PM
#Compare 4 and 8 tags experiments

library(ggplot2)
library(optimTE)


CRD232Tag4.des = optCRD(nTrt = 2, bRep  = 3, tRep  = 2, nPlot = 4, iter  = 1000)

CRD232Tag4 = simDataChisq(design =  CRD232Tag4.des,
                          blk.str2 = "Run", blk.str1 = "Ani",
                          trt.str = "Tag + Trt",
                          gamma.R = c(1, 5, 20),
                          gamma.A =c(1, 2, 10),
                          VC.resid = 1,
                          nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE)

CRD232Tag4New = simDataChisqNew(design =  CRD232Tag4.des, getVcEDF = getVcEDFEx1,
                          blk.str2 = "Run", blk.str1 = "Ani",
                          trt.str = "Tag + Trt",
                          gamma.R = c(1, 5, 20),
                          gamma.A =c(1, 2, 10),
                          VC.resid = 1,
                          nS = 4, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE)


(aov.table =summaryAovTwoPhase(test,  blk.str2 = "Run", blk.str1 = "Ani",
                              trt.str = "Tag + Trt")$A)
(rowNames = extractName(rownames(aov.table)))


CRD232Tag4New.EDF = simDataChisqNew(design =  CRD232Tag4.des, getVcEDF = getVcEDFEx1,
                              blk.str2 = "Run", blk.str1 = "Ani",
                              trt.str = "Tag + Trt",
                              gamma.R = c(0, 0.25, 1, 4, 100),
                              gamma.A = c(10^((-8:8)/2)),
                              VC.resid = 1,
                              nS = 4, nVc = 3, nSim = 100, seed = 1, neg.VC = TRUE)
edfPlot(CRD232Tag4New.EDF)


CRD232Tag4.neg = simDataChisq(design =  CRD232Tag4.des,
                          blk.str2 = "Run", blk.str1 = "Ani",
                          trt.str = "Tag + Trt",
                          gamma.R = c(1, 5, 20),
                          gamma.A =c(1, 2, 10),
                          VC.resid = 1,
                          nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = FALSE)

test <- data.frame(Case = rep(1:9, each = 3), CRD232Tag4$tempV[,c(5,4,3)])
test$VC <- factor(rownames(test)[1:3], rownames(test)[1:3])




xtable(dcast( melt(test, id.vars = c("Case", "VC")) , Case ~ variable +VC), digits = 4)


CRD232Tag4$tempV
CRD232Tag4.neg$tempV

CRD232Tag4$tempS
CRD232Tag4.neg$tempS




save.image("working07052017.Rdata")

CRD232Tag4.low = simDataChisq(design =  CRD232Tag4.des,
                          blk.str2 = "Run", blk.str1 = "Ani",
                          trt.str = "Tag + Trt",
                          gamma.R = c(1, 5, 20),
                          gamma.A =c(1, 2, 10),
                          VC.resid = 0.1,
                          nS = 4, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE)

test <- CRD232Tag4.low$tempV

c(test$REML.VC - test$REAL.VC)
c(test$LC.VC - test$REAL.VC)


CRD232Tag4.VC = simDataChisq(design =  CRD232Tag4.des,
                              blk.str2 = "Run", blk.str1 = "Ani",
                              trt.str = "Tag + Trt",
                              gamma.R = c(0, 1, 100),
                              gamma.A = c(10^((-1:1)*4)),
                              VC.resid = 1,
                              nS = 4, nVc = 3, nSim = 1000, seed = 12, neg.VC = TRUE)


test <- data.frame(Case = rep(1:9, each = 3), CRD232Tag4.VC$tempV[,c(5,4,3)])
test$VC <- factor(rownames(test)[1:3], rownames(test)[1:3])

xtable(dcast( melt(test, id.vars = c("Case", "VC")) , Case ~ variable +VC), digits = 4)


old = proc.time()
CRD232Tag4.EDF = simDataChisq(design =  CRD232Tag4.des,
                              blk.str2 = "Run", blk.str1 = "Ani",
                              trt.str = "Tag + Trt",
                              gamma.R = c(0, 0.25, 1, 4, 100),
                              gamma.A = c(10^((-8:8)/2)),
                              VC.resid = 1,
                              nS = 4, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE)
proc.time() - old

aov.table =summaryAovTwoPhase(CRD232Tag4.des,  blk.str2 = "Run", blk.str1 = "Ani",
                              trt.str = "Tag + Trt")$A
(rowNames = extractName(rownames(aov.table)))


xx =  simDataChisq(design =  CRD232Tag4.des,
                   blk.str2 = "Run", blk.str1 = "Ani",
                   trt.str = "Tag + Trt",
                   gamma.R = NA,
                   gamma.A =c(10^((-8:8)/2)),
                   VC.resid = 1,
                   nS = 2, nVc = 2, nSim = 10, seed = 1, 
                   row.MS = rowNames[c(8, 11)])
CRD232Tag4.EDF$tempS = rbind(CRD232Tag4.EDF$tempS, xx$tempS)


CRD232Tag4.EDF$tempS<- CRD232Tag4.EDF$tempS[,-5] 

edfPlot(CRD232Tag4.EDF)



old = proc.time()
CRD262Tag4.des = optCRD(nTrt = 2, bRep  = 6, tRep  = 2, nPlot = 4, iter  = 1000)

design.summary.CRD(CRD262Tag4.des, FALSE)


CRD262Tag4 = simDataChisq(design =  CRD262Tag4.des,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A = c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE)

aov.table =summaryAovTwoPhase(CRD262Tag4.des,  blk.str2 = "Run", blk.str1 = "Ani",
                              trt.str = "Tag + Trt")$A
(rowNames = extractName(rownames(aov.table)))


xx =  simDataChisq(design =  CRD262Tag4.des,
                    blk.str2 = "Run", blk.str1 = "Ani",
                    trt.str = "Tag + Trt",
                    gamma.R = NA,
                    gamma.A =c(10^((-8:8)/2)),
                    VC.resid = 1,
                    nS = 2, nVc = 2, nSim = 100, seed = 1, 
                    row.MS = rowNames[c(8, 11)])
CRD262Tag4$tempS = rbind(CRD262Tag4$tempS, xx$tempS)

# pdf(width = 12, file = "CRD232Tag4.pdf")
# edfPlot(CRD262Tag4)
# dev.off()

 old = proc.time()
CRD262Tag8.des = optCRD(nTrt = 2, bRep  = 6, tRep  = 2, nPlot = 8, iter  = 1000)
proc.time() - old

design.summary.CRD(design1.df, FALSE)

CRD262Tag8 = simDataChisq(design =  CRD262Tag8.des,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE)


aov.table =summaryAovTwoPhase(CRD262Tag8.des,  blk.str2 = "Run", blk.str1 = "Ani",
                              trt.str = "Tag + Trt")$A
(rowNames = extractName(rownames(aov.table)))


xx =  simDataChisq(design =  CRD262Tag8.des,
                   blk.str2 = "Run", blk.str1 = "Ani",
                   trt.str = "Tag + Trt",
                   gamma.R = NA,
                   gamma.A =c(10^((-8:8)/2)),
                   VC.resid = 1,
                   nS = 2, nVc = 2, nSim = 100, seed = 1, 
                   row.MS = rowNames[c(8, 11)])
CRD262Tag8$tempS = rbind(CRD262Tag8$tempS, xx$tempS)

proc.time() - old

#edfPlot(sim8Tags)

pdf(width = 12, file = "CRD262Tag4vsTag8.pdf")
edfPlotCompare(sim4Tags, CRD262Tag8,  compare = c("4Tag", "8Tag"))
dev.off()

edfPlotCompare.REML(CRD262Tag4, CRD262Tag8,  compare = c("4Tag", "8Tag"))


edfPlotCompare(sim4Tags, sim8Tags,  compare = c("4Tag", "8Tag"), method = "REAL")

save.image("Working.Rdata")


###############################################################################
#Compare 4 and 8 tags experiments    CRD 2

old = proc.time()
CRD462Tag4.des = optCRD(nTrt = 4, bRep  = 6, tRep  = 2, nPlot = 4, iter  = 1000)
proc.time() - old

design.summary.CRD(CRD462Tag4.des, FALSE)


old = proc.time()  
CRD462Tag4 = simDataChisq(design =  CRD462Tag4.des,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE)

(aov.table = summaryAovTwoPhase(CRD462Tag4.des,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(CRD462Tag4.des)))$A)

(rowNames = extractName(rownames(aov.table)))

xx =  simDataChisq(design =  CRD462Tag4.des,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10, seed = 1, 
            row.MS = rowNames[c(8, 11)])

CRD462Tag4$tempS = rbind(CRD462Tag4$tempS, xx$tempS)


old = proc.time()
CRD462Tag8.des = optCRD(nTrt = 4, bRep  = 6, tRep  = 2,  nPlot = 8, iter  = 1000)
proc.time() - old

design.summary.CRD(CRD462Tag8.des, FALSE)

aov.table =summaryAovTwoPhase(CRD462Tag8.des,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt")$A

rowNames = extractName(rownames(aov.table))

CRD462Tag8 = simDataChisq(design =  CRD462Tag8.des,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 1000, seed = 10, neg.VC = TRUE)


(aov.table = summaryAovTwoPhase(CRD462Tag8.des,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(CRD462Tag8.des)))$A)

(rowNames = extractName(rownames(aov.table)))

xx =  simDataChisq(design =  CRD462Tag8.des,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10, seed = 1, 
            row.MS = rowNames[c(8, 11)])

CRD462Tag8$tempS = rbind(CRD462Tag8$tempS, xx$tempS)

edfPlotCompare.REML(CRD462Tag4, CRD462Tag8,  compare = c("4Tag", "8Tag"))

edfPlot(CRD462Tag8)

pdf(width = 12, file = "CRD462Tag4vsTag8.pdf")
edfPlotCompare(sim4Tags, sim8Tags,  compare = c("4Tag", "8Tag"))
dev.off()

edfPlotCompare(CRD462Tag4, CRD462Tag8,  compare = c("4Tag", "8Tag"), method = "REAL")

save.image("CRD462Tag4vsTag8.Rdata")

summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", latex = TRUE, fixed.names = c("\\gamma", "\\tau"))


summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", latex = TRUE, fixed.names = c("\\gamma", "\\tau"))


###############################################################################
#Compare 4 and 8 tags experiments    CRD 3

old = proc.time()
CRD822Tag4.des = optCRD(nTrt = 8, bRep  = 2, tRep  = 2, nPlot = 4, iter  = 1000)
proc.time() - old

design.summary.CRD(CRD822Tag4.des, FALSE)

(aov.table = summaryAovTwoPhase(CRD822Tag4.des,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(CRD822Tag4.des)))$A)

(rowNames = extractName(rownames(aov.table)))

old = proc.time()  
CRD822Tag4 = simDataChisq(design =  CRD822Tag4.des,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE, checkFun = FALSE)
proc.time() - old 


(aov.table = summaryAovTwoPhase(CRD822Tag4.des,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(CRD822Tag4.des)))$A)

(rowNames = extractName(rownames(aov.table)))

xx =  simDataChisq(design =  CRD822Tag4.des,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10, seed = 1, 
            row.MS = rowNames[c(9, 12)])

CRD822Tag4$tempS = rbind(CRD822Tag4$tempS, xx$tempS)

            
edfPlot(CRD822Tag4)

old = proc.time()
CRD822Tag8.des = optCRD(nTrt = 8, bRep  = 2, tRep  = 2,  nPlot = 8, iter  = 1000)
proc.time() - old

design.summary.CRD(design1.df, FALSE)

aov.table = summaryAovTwoPhase(CRD822Tag8.des,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt")$A

(rowNames = extractName(rownames(aov.table)))

CRD822Tag8 = simDataChisq(design =  CRD822Tag8.des,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE)


(aov.table = summaryAovTwoPhase(CRD822Tag8.des,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(CRD822Tag8.des)))$A)

(rowNames = extractName(rownames(aov.table)))

xx =  simDataChisq(design =  CRD822Tag8.des,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10, seed = 1, 
            row.MS = rowNames[c(8, 11)])

CRD822Tag8$tempS = rbind(CRD822Tag8$tempS, xx$tempS)


edfPlot(sim8Tags)

pdf(width = 12, file = "CRD822Tag4vsTag8.pdf")
edfPlotCompare.REML(CRD822Tag4, CRD822Tag8,  compare = c("4Tag", "8Tag"))
dev.off()

edfPlotCompare(sim4Tags, sim8Tags,  compare = c("4Tag", "8Tag"), method = "REAL")

save.image("CRD822Tag4vsTag8.Rdata")

###############################################################################
#Compare 4 and 8 tags experiments RBD

design.df = optRBD(nTrt = 4, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 4, iter = 1000,
 confoundCag = FALSE, upperValue = 1)


 design1.df = optRBD(nTrt = 4, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 4, iter = 1000,
 confoundCag = TRUE, upperValue = 1)

  design.summary.RBD(design.df, FALSE)
 design.summary.RBD(design1.df, FALSE)

 design2.df = optRBD(nTrt = 4, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 8, iter = 100,
 confoundCag = FALSE, upperValue = 1)

 design3.df = optRBD(nTrt = 4, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 8, iter = 1000,
 confoundCag = TRUE, upperValue = 1)

   design.summary.RBD(design2.df, FALSE)
 design.summary.RBD(design3.df, FALSE)

 summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A

summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A

summaryAovTwoPhase(design2.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A

summaryAovTwoPhase(design3.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A

###############################################################################
#Compare 4 and 8 tags experiments RBD   1
library(snowfall)

 old = proc.time()
design.df = optRBD(nTrt = 4, bRep  = 4, nCag  = 4, tRep  = 2, nPlot = 4, iter = 1000,
 confoundCag = TRUE, upperValue = 1)

proc.time() - old

design.summary.RBD(design.df, FALSE)

(aov.table = summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

( rowNames = extractName(rownames(aov.table)))


sim4Tags1 = simDataChisqPRBD(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 4, nSim = 10, seed = 1, neg.VC = TRUE, checkFun = FALSE)

edfPlot1.comCag(sim4Tags1)

sim4Tags0 = sim4Tags

sim4Tags0$tempS = sim4Tags0$tempS[which(sim4Tags0$tempS$gamma.cag == 0.25), -2]

edfPlot(sim4Tags0, "Within RunBetweenCag(Ani)Residual")

sim4Tags1 = sim4Tags

sim4Tags1$tempS = sim4Tags1$tempS[which(sim4Tags1$tempS$gamma.cag == 100), -2]

edfPlot(sim4Tags1, "Within RunBetweenCag(Ani)Residual")

edfPlotCompare(sim4Tags0, sim4Tags1, compare = c("0.25", "100"),
method = "REML", rowNames = "Within RunBetweenCag:AniResidual")

#Run is fixed effect

(aov.table =summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

( rowNames = extractName(rownames(aov.table)))

sim4Tags1RunCagFix = simDataChisqPRBD(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(8, 11)])

edfPlot1.comCag(sim4Tags1RunCagFix)

sim4Tags1CagFix = simDataChisqPRBD(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0,  0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 3, nSim = 10, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(3, 8, 11)])

sim4Tags1.all = sim4Tags1

sim4Tags1.all$tempS = rbind( sim4Tags1$tempS, 
                             sim4Tags1CagFix$tempS, 
                             sim4Tags1RunCagFix$tempS)




 (aov.table = summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

( rowNames = extractName(rownames(aov.table)))


sim4Tags2 = simDataChisqPRBD(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 6, nVc = 4, nSim = 10, seed = 1, neg.VC = TRUE, checkFun = FALSE)

 (aov.table = summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

( rowNames = extractName(rownames(aov.table)))

sim4Tags2RunFix = simDataChisqPRBD(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(8, 11, 14)])

sim4Tags2CagFix = simDataChisqPRBD(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0,  0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(3, 4, 11, 14)])

sim4Tags2CagRunFix = simDataChisqPRBD(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(11, 14)])

sim4Tags2.all = sim4Tags2

sim4Tags2.all$tempS = rbind( sim4Tags2$tempS, 
                             sim4Tags2RunFix$tempS, 
                             sim4Tags2CagFix$tempS, 
                             sim4Tags2CagRunFix$tempS)

edfPlot1.comCag(sim4Tags2.all)  + geom_abline(intercept = 8, slope = 0)

sim4Tags2.comp = sim4Tags2.all
 sim4Tags1.comp =  sim4Tags1.all
 
 sim4Tags1.comp$tempS =  sim4Tags1.comp$tempS[which(is.na(sim4Tags1.comp$tempS$gamma.cag)), ]

 sim4Tags1.comp$tempS$gamma.cag = as.factor("CagConfRun")


 sim4Tags2.comp$tempS = rbind( sim4Tags2.comp$tempS, 
                             sim4Tags1.comp$tempS)

pdf("RBD442Tag4EDF.pdf", width = 12)
  edfPlot2.comCag(sim4Tags2.comp)
dev.off()

 sim4Tags2.run0  = sim4Tags2.comp
 
 sim4Tags2.run0$tempS =  sim4Tags2.run0$tempS[which(sim4Tags2.run0$tempS$gamma.run == 0), ]

pdf("RBD442Tag4EDFRun0.pdf", width = 7)
edfPlot2.comCag(sim4Tags2.run0)
dev.off()



negVCsim4Tags2 = simDataChisqPRBD(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 6, nVc = 4, nSim = 100, seed = 1, neg.VC = FALSE, checkFun = FALSE)

edfPlot1.comCag(negVCsim4Tags2)

 (aov.table = summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

( rowNames = extractName(rownames(aov.table)))

negVCsim4Tags2RunFix = simDataChisqPRBD(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 10000, seed = 1, neg.VC = FALSE,
                 row.MS =rowNames[c(8, 11, 14)])

negVCsim4Tags2CagFix = simDataChisqPRBD(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0,  0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = FALSE,
                 row.MS =rowNames[c(3, 4, 11, 14)])

negVCsim4Tags2.all = negVCsim4Tags2

negVCsim4Tags2.all$tempS = rbind( negVCsim4Tags2$tempS, 
                             negVCsim4Tags2RunFix$tempS, 
                             negVCsim4Tags2CagFix$tempS)
                             
 (aov.table = summaryAovTwoPhase(design2.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

( rowNames = extractName(rownames(aov.table)))

sim8Tags1 = simDataChisqPRBD(design =  design2.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 5, nVc = 4, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)


sim8Tags1RunFix  = simDataChisqPRBD(design =  design2.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 100, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(3, 7, 11, 14)])

sim8Tags1CagFix  = simDataChisqPRBD(design =  design2.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R =  c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(3, 11, 14)])

sim8Tags1RunCagFix  = simDataChisqPRBD(design =  design2.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R =  NA,
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 100, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(11, 14)])

sim8Tags1.all = sim8Tags1

sim8Tags1.all$tempS = rbind( sim8Tags1$tempS, 
                             sim8Tags1RunFix$tempS, 
                             sim8Tags1CagFix$tempS, 
                             sim8Tags1RunCagFix$tempS)

 (aov.table = summaryAovTwoPhase(design2.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

( rowNames = extractName(rownames(aov.table)))

sim8Tags1.Cag = simDataChisqPRBDCag(design =  design2.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 5, nVc = 4, nSim = 1000, seed = 1, neg.VC = TRUE, checkFun = FALSE)

sim8Tags1.Cag.fixRun = simDataChisqPRBDCag(design =  design2.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 100, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(7, 11, 14)])


sim8Tags1.Cag.all = sim8Tags1.Cag

sim8Tags1.Cag.all$tempS = rbind( sim8Tags1.Cag$tempS, 
                             sim8Tags1.Cag.fixRun$tempS)

                 
                 
# pdf("RBD442Tag8CagEDF.pdf", width = 12)
# edfPlot1.comCag(sim8Tags1.Cag.all, rowNames = ".")
# dev.off()

 (aov.table = summaryAovTwoPhase(design3.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

( rowNames = extractName(rownames(aov.table)))

 sim8Tags2 = simDataChisqPRBD(design =  design3.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)

sim8Tags2RunCagFix  = simDataChisqPRBD(design =  design3.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R =NA,
            gamma.C =NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(9, 12)])


sim8Tags2.all = sim8Tags2

sim8Tags2.all$tempS = rbind( sim8Tags2$tempS, 
                             sim8Tags2RunCagFix$tempS)

sim8Tags2.all$tempS = sim8Tags2.all$tempS[,-2]

 sim4Tags2.cag0  = sim4Tags2.comp
 sim4Tags2.cag0$tempS =  sim4Tags2.cag0$tempS[which(sim4Tags2.cag0$tempS$gamma.cag == 0), ]
sim4Tags2.cag0$tempS = sim4Tags2.cag0$tempS[,-2]


sim4Tags1.cag0  = sim4Tags1.all
sim4Tags1.cag0$tempS =  sim4Tags1.cag0$tempS[which(sim4Tags1.cag0$tempS$gamma.cag == 0), ]
sim4Tags1.cag0$tempS = sim4Tags1.cag0$tempS[,-2]

xx = data.frame(as.factor(NA), as.factor(c(10^((-8:8)/2))), 
          matrix(8, nrow = length(c(10^((-8:8)/2))), ncol = 3))
          
colnames(xx) = 
colnames(sim4Tags1.cag0$tempS)

sim4Tags1.cag0$tempS = rbind(sim4Tags1.cag0$tempS, xx)


sim8Tags1.cag0  = sim8Tags1.all
sim8Tags1.cag0$tempS =  sim8Tags1.cag0$tempS[which(sim8Tags1.cag0$tempS$gamma.cag == 0), ]
sim8Tags1.cag0$tempS = sim8Tags1.cag0$tempS[,-2]

xx = data.frame(as.factor(NA), as.factor(c(10^((-8:8)/2))), 
          matrix(7, nrow = length(c(10^((-8:8)/2))), ncol = 3))
          
colnames(xx) = 
colnames(sim8Tags1.cag0$tempS)

sim8Tags1.cag0$tempS = rbind(sim8Tags1.cag0$tempS, xx)


pdf("RBD442Tag4vsTag8.pdf", width = 12)
edfPlotCompare1(sim4Tags1.cag0, sim4Tags2.cag0, sim8Tags1.cag0, sim8Tags2.all,  
compare = c("4TagCagConfRun", "4TagCagConfTag", "8TagCagConfRun", "8TagCagConfTag"), 
rowNames = ".") 
dev.off()





save.image("RBD442Tag4vsTag8.Rdata")

###############################################################################
#Compare 4 and 8 tags experiments BIBD




BIBDdesign4 = optBIBD(nTrt = 4, bRep  = 3, nCag  = 4, tRep  = 2, nPlot = 4, iter = 1000,
    confoundCag = TRUE)


design.summary.RBD(BIBDdesign4)

(aov.table =summaryAovTwoPhase(BIBDdesign4,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))

sim4Tags = simDataChisq1(design = BIBDdesign4,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE, checkFun = FALSE)


xx = simDataChisq1(design = BIBDdesign4,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 1000,  seed = 1,
            row.MS =rowNames[c(10, 13)])


sim4Tags1$tempS = rbind(sim4Tags$tempS, xx$tempS)
sim4Tags1$tempS = sim4Tags1$tempS[,-2]


edfPlot(sim4Tags1, "Within RunBetweenCag:AniResidual")



BIBDdesign48 = optBIBD(nTrt = 4, bRep  = 3, nCag  = 4, tRep  = 2, nPlot = 8, iter = 1000,
    confoundCag = TRUE)

design.summary.RBD(BIBDdesign48)

(aov.table = summaryAovTwoPhase(BIBDdesign48,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))

sim8Tags = simDataChisq1(design = BIBDdesign48,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE, checkFun = FALSE)

xx = simDataChisq1(design = BIBDdesign48,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 1000,  seed = 1,
            row.MS =rowNames[c(9, 12)])

sim8Tags1$tempS = rbind(sim8Tags$tempS, xx$tempS)
sim8Tags1$tempS = sim8Tags1$tempS[,-2]


edfPlot(sim8Tags1, "Within RunBetweenCag:AniResidual")


 
pdf(width = 12, file = "BIBD432Tag4vsTag8.pdf")
edfPlotCompare(sim4Tags2.cag0, sim8Tags1,  compare = c("4Tag", "8Tag"), rowNames = ".")
dev.off()



#############################################################################

BIBDdesign4 = optBIBD(nTrt = 5, bRep  = 4, nCag  = 5, tRep  = 2, nPlot = 4, iter = 1000,
    confoundCag = FALSE)

 design.summary.RBD(BIBDdesign4, FALSE)

 BIBDdesign4$ani = interaction(BIBDdesign4$Cag, BIBDdesign4$Ani)


  gamma.run = 10
    gamma.cag = 5
  gamma.ani = 2
  VC.resid = 1
  

     run.eff = rnorm(nlevels(BIBDdesign4$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
    cag.eff = rnorm(nlevels(BIBDdesign4$Cag), mean = 0, sd = sqrt(gamma.cag * VC.resid))
    ani.eff = rnorm(nlevels(BIBDdesign4$ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
    trt.eff = runif(nlevels(BIBDdesign4$Trt), 0, 10)
    tag.eff = runif(nlevels(BIBDdesign4$Tag), 0, 1)
    res.eff = rnorm(nrow(BIBDdesign4), mean = 0, sd = sqrt(VC.resid))
    gm = 10

      y = gm + with(BIBDdesign4, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] +
            tag.eff[Tag] + trt.eff[Trt]) + res.eff

     summaryAovTwoPhase(BIBDdesign4,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt", response = y)


   Tag.mat = with(BIBDdesign4, as.matrix(table(1:nrow(BIBDdesign4), Tag)))
Trt.mat = with(BIBDdesign4, as.matrix(table(1:nrow(BIBDdesign4), Trt)))
   
X = cbind(1, Trt.mat, Tag.mat)


     summaryAovTwoPhase(BIBDdesign4,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt", response = y -  (trtProj(BIBDdesign4) %*% y)[BIBDdesign4$Trt])

