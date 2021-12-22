#22/05/2013 02:26:39 PM
#Compare 4 and 8 tags experiments



 old = proc.time()
design.df = optCRD(nTrt = 2, bRep  = 4, tRep  = 2, nPlot = 4, iter  = 1000)
proc.time() - old

design.summary.CRD(design.df, FALSE)

sim4Tags = simDataChisq(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE)

pdf(width = 12, file = "CRD232Tag4.pdf")
edfPlot(sim4Tags)
dev.off()

 old = proc.time()
design1.df = optCRD(nTrt = 2, bRep  = 6, tRep  = 2, nPlot = 8, iter  = 1000)
proc.time() - old

design.summary.CRD(design1.df, FALSE)

sim8Tags = simDataChisq(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE)


edfPlot(sim8Tags)

pdf(width = 12, file = "CRD262Tag4vsTag8.pdf")
edfPlotCompare(sim4Tags, sim8Tags,  compare = c("4Tag", "8Tag"))
dev.off()

edfPlotCompare(sim4Tags, sim8Tags,  compare = c("4Tag", "8Tag"), method = "REAL")


###############################################################################
#Compare 4 and 8 tags experiments    CRD



 old = proc.time()
design.df = optCRD(nTrt = 3, bRep  = 12, tRep  = 2, nPlot = 4, iter  = 1000)
proc.time() - old

design.summary.CRD(design.df, FALSE)

aov.table = summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(design.df)))

summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", latex = TRUE, fixed.names = c("\\gamma", "\\tau"))

rowNames = extractName(rownames(aov.table))

sim4Tags = simDataChisq(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE, checkFun = FALSE)
            
edfPlot(sim4Tags)

old = proc.time()
design1.df = optCRD(nTrt = 3, bRep  = 12, tRep  = 2, nPlot = 8, iter  = 1000)
proc.time() - old

design.summary.CRD(design1.df, FALSE)

aov.table =summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt")$A

rowNames = extractName(rownames(aov.table))

 summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", latex = TRUE, fixed.names = c("\\gamma", "\\tau"))

sim8Tags = simDataChisq(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE)


            row.MS = rowNames[c(3,4, 9, 12)]



edfPlot(sim8Tags)

edfPlotCompare(sim4Tags, sim8Tags,  compare = c("4Tag", "8Tag"), method = "REML")

edfPlotCompare(sim4Tags, sim8Tags,  compare = c("4Tag", "8Tag"), method = "REAL")

###############################################################################
#Compare 4 and 8 tags experiments    CRD 2

old = proc.time()
design.df = optCRD(nTrt = 4, bRep  = 6, tRep  = 2, nPlot = 4, iter  = 1000)
proc.time() - old

design.summary.CRD(design.df, FALSE)

aov.table = summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(design.df)))

rowNames = extractName(rownames(aov.table))

old = proc.time()  
sim4Tags = simDataChisqP(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)
proc.time() - old 
  sfStop()


(aov.table = summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(design.df)))$A)

(rowNames = extractName(rownames(aov.table)))

xx =  simDataChisqP(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10000, seed = 1, 
            row.MS = rowNames[c(8, 11)])

sim4Tags$tempS = rbind(sim4Tags$tempS, xx$tempS)

pdf(width = 12, file = "CRD462Tag4.pdf")
edfPlot(sim4Tags)
dev.off()
            
edfPlot(sim4Tags)

old = proc.time()
design1.df = optCRD(nTrt = 4, bRep  = 6, tRep  = 2,  nPlot = 8, iter  = 1000)
proc.time() - old

design.summary.CRD(design1.df, FALSE)

aov.table =summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt")$A

rowNames = extractName(rownames(aov.table))

sim8Tags = simDataChisqP(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE)


(aov.table = summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(design.df)))$A)

(rowNames = extractName(rownames(aov.table)))

xx =  simDataChisqP(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10000, seed = 1, 
            row.MS = rowNames[c(8, 11)])

sim8Tags$tempS = rbind(sim8Tags$tempS, xx$tempS)


edfPlot(sim8Tags)

pdf(width = 12, file = "CRD462Tag4vsTag8.pdf")
edfPlotCompare(sim4Tags, sim8Tags,  compare = c("4Tag", "8Tag"))
dev.off()

edfPlotCompare(sim4Tags, sim8Tags,  compare = c("4Tag", "8Tag"), method = "REAL")

save.image("CRD462Tag4vsTag8.Rdata")

summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", latex = TRUE, fixed.names = c("\\gamma", "\\tau"))


summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", latex = TRUE, fixed.names = c("\\gamma", "\\tau"))


###############################################################################
#Compare 4 and 8 tags experiments    CRD 3

old = proc.time()
design.df = optCRD(nTrt = 8, bRep  = 2, tRep  = 2, nPlot = 4, iter  = 1000)
proc.time() - old

design.summary.CRD(design.df, FALSE)

(aov.table = summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(design.df)))$A)

(rowNames = extractName(rownames(aov.table)))

old = proc.time()  
sim4Tags = simDataChisqP(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)
proc.time() - old 
  sfStop()


(aov.table = summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(design.df)))$A)

(rowNames = extractName(rownames(aov.table)))

xx =  simDataChisqP(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10000, seed = 1, 
            row.MS = rowNames[c(9, 12)])

sim4Tags$tempS = rbind(sim4Tags$tempS, xx$tempS)

            
edfPlot(sim4Tags)

old = proc.time()
design1.df = optCRD(nTrt = 8, bRep  = 2, tRep  = 2,  nPlot = 8, iter  = 1000)
proc.time() - old

design.summary.CRD(design1.df, FALSE)

aov.table =summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt")$A

(rowNames = extractName(rownames(aov.table)))

sim8Tags = simDataChisqP(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE)


(aov.table = summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(design.df)))$A)

(rowNames = extractName(rownames(aov.table)))

xx =  simDataChisqP(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10000, seed = 1, 
            row.MS = rowNames[c(8, 11)])

sim8Tags$tempS = rbind(sim8Tags$tempS, xx$tempS)


edfPlot(sim8Tags)

pdf(width = 12, file = "CRD822Tag4vsTag8.pdf")
edfPlotCompare(sim4Tags, sim8Tags,  compare = c("4Tag", "8Tag"))
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
            nS = 4, nVc = 4, nSim = 1000, seed = 1, neg.VC = TRUE, checkFun = FALSE)

edfPlot1.comCag(sim4Tags)

sim4Tags0 = sim4Tags

sim4Tags0$tempS = sim4Tags0$tempS[which(sim4Tags0$tempS$gamma.cag == 0.25), -2]

edfPlot(sim4Tags0, "Within RunBetweenCag:AniResidual")

sim4Tags1 = sim4Tags

sim4Tags1$tempS = sim4Tags1$tempS[which(sim4Tags1$tempS$gamma.cag == 100), -2]

edfPlot(sim4Tags1, "Within RunBetweenCag:AniResidual")

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
            nS = 2, nVc = 2, nSim = 1000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(8, 11)])

edfPlot1.comCag(sim4Tags1RunCagFix)

sim4Tags1CagFix = simDataChisqPRBD(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0,  0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 3, nSim = 1000, seed = 1, neg.VC = TRUE,
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
            nS = 6, nVc = 4, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)

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

sim4Tags2.all = sim4Tags2

sim4Tags2.all$tempS = rbind( sim4Tags2$tempS, 
                             sim4Tags2RunFix$tempS, 
                             sim4Tags2CagFix$tempS)


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
            gamma.R =NA,
            gamma.C = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(7, 11, 14)])

sim8Tags1CagFix  = simDataChisqPRBD(design =  design2.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R =  c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(3, 11, 14)])


sim8Tags1.all = sim8Tags1

sim8Tags1.all$tempS = rbind( sim8Tags1$tempS, 
                             sim8Tags1RunFix$tempS, 
                             sim8Tags1CagFix$tempS)


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

save.image("RBD442Tag4vsTag8.Rdata")

edfPlotCompare(sim4Tags1, sim8Tags,  compare = c("4Tag", "8Tag"),
method = "REML", rowNames = "Within RunBetweenCag:AniResidual") +
 geom_abline (intercept = 8, slope=0)

summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt", latex = TRUE, fixed.names = c("\\gamma", "\\tau"))

summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt", latex = TRUE, fixed.names = c("\\gamma", "\\tau"))

summaryAovTwoPhase(design2.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt", latex = TRUE, fixed.names = c("\\gamma", "\\tau"))

summaryAovTwoPhase(design3.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt", latex = TRUE, fixed.names = c("\\gamma", "\\tau"))



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
edfPlotCompare(sim4Tags1, sim8Tags1,  compare = c("4Tag", "8Tag"),
method = "REML", rowNames = "Within RunBetweenCag:AniResidual")
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

