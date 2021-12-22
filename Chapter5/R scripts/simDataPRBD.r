

simDataChisqPRBDCag =
function(design, blk.str2, blk.str1, trt.str, gamma.R,  gamma.C, gamma.A, 
                  VC.resid, nS, nVc, nSim, seed, row.MS = NA, neg.VC = TRUE, checkFun = FALSE){
  
  sfInit(parallel = TRUE, cpus = 4)

  aov.table = summaryAovTwoPhase(design.df = design, blk.str1 =  blk.str1,
                          blk.str2 = blk.str2, trt.str = trt.str, response = rnorm(nrow(design)))
   
     tempS = data.frame()

 
# create progress bar
  for(i in 1:length(gamma.R)){
    gamma.run = gamma.R[i]

    cat("gamma.run = ", gamma.run, "\n")
    
    for(k in 1:length(gamma.C)){
      gamma.cag = gamma.C[k]

  cat("gamma.cag = ", gamma.cag, "\n")
    
  
     test = function(j) {
  
      gamma.ani = gamma.A[j]

      REML.EDF = numeric(nS); LC.EDF = numeric(nS); REAL.EDF  = numeric(nS);
 
      real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid), (gamma.run * VC.resid))  
      index =  1:length(real.VC)
      if(any(is.na(real.VC))){
         index = index[-which(is.na(real.VC))]
         real.VC =  real.VC[-which(is.na(real.VC))]
             
      }
  
      set.seed(seed)
      simMS = suppressWarnings(apply(aov.table$A, 1, 
          function(x) sum(fracToNum(x[index+1]) * real.VC) * rchisq(nSim, fracToNum(x[1]))/fracToNum(x[1])))
 
     
        counter = 0
       #simulation
       while( counter < nSim){
       
        tmp = try(getVcEDF(aov.table = aov.table, MS = simMS[counter+1, ], true.VC = real.VC, 
                    neg.VC = neg.VC, row.MS = row.MS), TRUE)
        

        if(class(tmp) =="try-error"){
          counter = counter + 1
          next
        }
        
        tmpS = tmp$S

        REML.EDF  =  rbind(REML.EDF,     tmpS[,"REML.EDF"] )
        LC.EDF    =  rbind(LC.EDF  ,     tmpS[,"LC.EDF"]   )
        REAL.EDF  =  rbind(REAL.EDF,     tmpS[,"TRUE.EDF"] )

        counter = counter + 1
      }
 
      
      REML.EDF = apply(REML.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
       LC.EDF = apply(LC.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
       REAL.EDF = apply(REAL.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
          
       return(cbind(gamma.run,  gamma.cag, gamma.ani, REML.EDF,  LC.EDF, REAL.EDF)["Within RunBetweenCagResidual",]) 
      
    }

    sfExport("nSim", "gamma.cag", "gamma.A", "VC.resid", "gamma.run", "nS", "fracToNum", "nVc", 
    "getVcEDF", "extractName", "row.MS", "neg.VC", "aov.table", "seed")
    
    sfLibrary(MASS)
    
    if(checkFun) browser()
 
     tempS = rbind(tempS, t(sfSapply(1:length(gamma.A), fun = test)))
     }
  }
   sfStop()

    tempS$gamma.run = as.factor(tempS$gamma.run)
    tempS$gamma.cag = as.factor(tempS$gamma.cag)
    tempS$gamma.ani  = as.factor(tempS$gamma.ani)
    rownames(tempS) = paste(rep("Within RunBetweenCagResidual", nrow(tempS)), 1:nrow(tempS), sep = "")
    return(list(tempS = tempS))     
  }

##########################################################################################

simDataChisqPRBD =
function(design, blk.str2, blk.str1, trt.str, gamma.R,  gamma.C, gamma.A, 
                  VC.resid, nS, nVc, nSim, seed, row.MS = NA, neg.VC = TRUE, checkFun = FALSE){
  
  #sfInit(parallel = TRUE, cpus = 4)

  aov.table = summaryAovTwoPhase(design.df = design, blk.str1 =  blk.str1,
                          blk.str2 = blk.str2, trt.str = trt.str, response = rnorm(nrow(design)))
   
     tempS = data.frame()

 
# create progress bar
  for(i in 1:length(gamma.R)){
    gamma.run = gamma.R[i]

    cat("gamma.run = ", gamma.run, "\n")
    
    for(k in 1:length(gamma.C)){
      gamma.cag = gamma.C[k]

  cat("gamma.cag = ", gamma.cag, "\n")
    
  
     test = 
       function(j) {
  
      gamma.ani = gamma.A[j]

      REML.EDF = numeric(nS); LC.EDF = numeric(nS); REAL.EDF  = numeric(nS);
 
      real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid), (gamma.run * VC.resid))  
      index =  1:length(real.VC)
      if(any(is.na(real.VC))){
         index = index[-which(is.na(real.VC))]
         real.VC =  real.VC[-which(is.na(real.VC))]
             
      }
  
      set.seed(seed)
      simMS = suppressWarnings(apply(aov.table$A, 1, 
          function(x) sum(fracToNum(x[index+1]) * real.VC) * rchisq(nSim, fracToNum(x[1]))/fracToNum(x[1])))
 
     
        counter = 0
       #simulation
       while( counter < nSim){
       
        tmp = try(getVcEDF(aov.table = aov.table, MS = simMS[counter+1, ], true.VC = real.VC, 
                    neg.VC = neg.VC, row.MS = row.MS), TRUE)
        

        if(class(tmp) =="try-error"){
          counter = counter + 1
          next
        }
        
        tmpS = tmp$S

        REML.EDF  =  rbind(REML.EDF,     tmpS[,"REML.EDF"] )
        LC.EDF    =  rbind(LC.EDF  ,     tmpS[,"LC.EDF"]   )
        REAL.EDF  =  rbind(REAL.EDF,     tmpS[,"TRUE.EDF"] )

        counter = counter + 1
      }
 
      
      REML.EDF = apply(REML.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
       LC.EDF = apply(LC.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
       REAL.EDF = apply(REAL.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
          
       return(cbind(gamma.run,  gamma.cag, gamma.ani, REML.EDF,  LC.EDF, REAL.EDF)["Within RunBetweenCag(Ani)Residual",]) 
      
    }

    #sfExport("nSim", "gamma.cag", "gamma.A", "VC.resid", "gamma.run", "nS", "fracToNum", "nVc", 
    #"getVcEDF", "extractName", "row.MS", "neg.VC", "aov.table", "seed")
    
    #sfLibrary(MASS)
    
    if(checkFun) browser()
 
     tempS = rbind(tempS, t(sapply(1:length(gamma.A), test)))
     }
  }
   #sfStop()

    tempS$gamma.run = as.factor(tempS$gamma.run)
    tempS$gamma.cag = as.factor(tempS$gamma.cag)
    tempS$gamma.ani  = as.factor(tempS$gamma.ani)
    rownames(tempS) = paste(rep("Within RunBetweenCag(Ani)Residual", nrow(tempS)), 1:nrow(tempS), sep = "")
    return(list(tempS = tempS))     
  }

##########################################################################################
# RBD cases  different initial designs 

design.df = optRBD(nTrt = 4, bRep  = 4, nCag  = 2, tRep  = 2, nPlot = 4, iter = 1000,
 confoundCag = TRUE, upperValue = 1)
 

design1.df = optRBD(nTrt = 4, bRep  = 4, nCag = 2, tRep  = 2, nPlot = 4, iter = 1000,
 confoundCag = FALSE, upperValue = 1)
  

 design.summary.RBD(design1.df, FALSE)


(aov.table =summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))

sim4Tags.ConRun = simDataChisqPRBD(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)


 sfStop()
 
 sim4Tags.ConRun$tempS = sim4Tags.ConRun$tempS[,-2]

edfPlot(sim4Tags.ConRun, "Within RunBetweenCag:AniResidual")

sim4Tags.ConRun.RunFix = simDataChisqPRBD(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(8, 11, 14)])
 sfStop()

edfPlot1(sim4Tags.ConRun.RunFix)



design1.df = optRBD(nTrt = 4, bRep  = 4, nCag  = 2, tRep  = 2, nPlot = 4, iter = 1000,
 confoundCag = FALSE, upperValue = 1)

 design.summary.RBD(design1.df, FALSE)

 (aov.table =summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))

sim4Tags.ConTag= simDataChisqPRBD(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 5, nVc = 4, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)
 sfStop()
 plots = edfPlot1(sim4Tags.ConTag)


sim4Tags.ConTag.CagFix= simDataChisqPRBD(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0,  0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE,
                 row.MS =rowNames[c(3,4, 9, 12)])


design2.df = optRBD(nTrt = 4, bRep  = 4, nCag  = 2, tRep  = 2, nPlot = 8, iter = 1000,
 confoundCag = TRUE, upperValue = 1)

 design.summary.RBD(design2.df, FALSE)

 (aov.table =summaryAovTwoPhase(design2.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))

sim8Tags.ConRun= simDataChisqPRBD(design =  design2.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)

 sfStop()
 

 design3.df = optRBD(nTrt = 4, bRep  = 4, nCag  = 2, tRep  = 2, nPlot = 8, iter = 1000,
 confoundCag = FALSE, upperValue = 1)
 
 design.summary.RBD(design3.df, FALSE)

  (aov.table =summaryAovTwoPhase(design3.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))


sim8Tags.ConTag= simDataChisqPRBD(design =  design3.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 4, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)

 plots = edfPlot1(sim8Tags.ConTag)

sim8Tags.ConTag.RunFix = simDataChisqPRBD(design =  design3.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(7, 11, 14)])
 sfStop()

edfPlot1(sim4Tags.ConRun.RunFix)


sim8Tags.ConTag.CagFix = simDataChisqPRBD(design =  design3.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(3, 11, 14)])

 sfStop()

sim8Tags.ConTag.CagFix$tempS = sim8Tags.ConTag.CagFix$tempS[,-2]

##########################################################################################
# RBD cases  different initial designs 

design.df = optRBD(nTrt =  2, bRep  = 8, nCag  = 4, tRep  = 2, nPlot = 4 iter = 1000,
 confoundCag = FALSE, upperValue = 1)
 


 design.summary.RBD(design.df, FALSE)


(aov.table =summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))

sim4Tags.ConRun = simDataChisqPRBD(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 6, nVc = 4, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)
 sfStop()
 
edfPlot1.comCag(sim4Tags.ConRun)

sim4Tags.ConRun.RunFix = simDataChisqPRBD(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(8, 11, 14)])
 sfStop()

edfPlot1.comCag(sim4Tags.ConRun.RunFix)


sim4Tags.ConRun.CagFix = simDataChisqPRBD(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(3, 4, 11, 14)])

 sfStop()

sim4Tags.ConRun.CagFix$tempS = sim4Tags.ConRun.CagFix$tempS[,-2]


edfPlot(sim4Tags.ConRun.CagFix, "Within RunBetweenCag:AniResidual")

sim4Tags.ConRun1 = sim4Tags.ConRun

sim4Tags.ConRun1$tempS = rbind( sim4Tags.ConRun$tempS, 
                                sim4Tags.ConRun.CagFix$tempS, 
                                sim4Tags.ConRun.RunFix$tempS)


pdf(width = 12, "RBD2842Tag4.pdf")
edfPlot1.comCag(sim4Tags.ConRun1) +
 geom_abline (intercept = 10, slope=0)
dev.off()

edfPlot1.comRun(sim4Tags.ConRun1)

design1.df = optRBD(nTrt = 2, bRep  = 8, nCag  = 4, tRep  = 2, nPlot = 4, iter = 1000,
 confoundCag = FALSE, upperValue = 1)

 design.summary.RBD(design1.df, FALSE)

 (aov.table =summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))

sim4Tags.ConTag= simDataChisqPRBD(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 4, nSim = 100, seed = 1, neg.VC = TRUE, checkFun = FALSE)
 sfStop()
edfPlot1.comCag(sim4Tags.ConTag)


design2.df = optRBD(nTrt = 2, bRep  = 8, nCag  = 4, tRep  = 2, nPlot = 8, iter = 1000,
 confoundCag = TRUE, upperValue = 1)

 design.summary.RBD(design2.df, FALSE)

 (aov.table = summaryAovTwoPhase(design2.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))

sim8Tags.ConRun = simDataChisqPRBD(design =  design2.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)
 sfStop()

 temp =  cbind(NA, NA, c(10^((-8:8)/2)), 10, 10, 10)
 
 colnames(temp) = colnames(sim8Tags.ConRun$tempS)
 
 temp1 = sim8Tags.ConRun
 temp1$tempS = rbind(sim8Tags.ConRun$tempS, temp)  
  

 design3.df = optRBD(nTrt = 2, bRep  = 8, nCag  = 4, tRep  = 2, nPlot = 8, iter = 1000,
 confoundCag = FALSE, upperValue = 1)
 
 design.summary.RBD(design3.df, FALSE)

  (aov.table =summaryAovTwoPhase(design3.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))


sim8Tags.ConTag= simDataChisqPRBD(design =  design3.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 5, nVc = 4, nSim = 100, seed = 1, neg.VC = TRUE, checkFun = FALSE)


sim8Tags.ConTag.RunFix = simDataChisqPRBD(design =  design3.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 100, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(7, 11, 14)])
 sfStop()


sim8Tags.ConTag1 = sim8Tags.ConTag

sim8Tags.ConTag1$tempS = rbind( sim8Tags.ConTag$tempS,
                                sim8Tags.ConTag.RunFix$tempS)

pdf(width = 12, "RBD2842TagTag8CagEDF.pdf")
edfPlot1.comCag(sim8Tags.ConTag1)
dev.off()

sim8Tags.ConTag.CagFix = simDataChisqPRBD(design =  design3.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(3, 11, 14)])

 sfStop()

sim8Tags.ConTag.CagFix$tempS = sim8Tags.ConTag.CagFix$tempS[,-2]


edfPlot1.comCag(temp1, ".")
temp1$tempS = temp1$tempS[,-2]

#comparing between the two designs 

temp2 = sim4Tags.ConRun1

temp2$tempS = temp2$tempS[which(temp2$tempS$gamma.cag == 0),-2]

pdf(width = 12, "RBD2842Tag4vsTag8.pdf")
edfPlotCompare(temp2, temp1, rowNames = ".", compare = c("4Tag", "8Tag"), method = "REML")
dev.off()

temp3 = sim4Tags.ConRun1
temp3$tempS = temp3$tempS[which(temp3$tempS$gamma.run == 0),]
edfPlot1.comCag(temp3, ".")

pdf("RBD2842Tag4run0.pdf")
edfPlot1.comCag(temp3, ".")
dev.off()


##########################################################################################
# RBD cases 

design.df = optRBD(nTrt = 2, bRep  = 8, nCag  = 4, tRep  = 2, nPlot = 4, iter = 1000,
 confoundCag = TRUE, upperValue = 1)

design.summary.RBD(design.df, FALSE)

design1.df = optRBD(nTrt = 2, bRep  = 8, nCag  = 6, tRep  = 2, nPlot = 4, iter = 1000,
 confoundCag = FALSE, upperValue = 1)

 design.summary.RBD(design1.df, FALSE)

(aov.table =summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))

sim4Tags = simDataChisqPRBD(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 6, nVc = 4, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)

 plots = edfPlot1(sim4Tags)

sim4TagsRunFix = simDataChisqPRBD(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = c(0,  0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 3, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(8, 11, 14)])


edfPlot1(sim4TagsRunFix)


sim4TagsCagFix = simDataChisqPRBD(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE,
                 row.MS =rowNames[c(3, 4, 11, 14)])

sim4TagsCagFix$tempS = sim4TagsCagFix$tempS[,-2]


edfPlot(sim4TagsCagFix, "Within RunBetweenCag:AniResidual")



design.summary.RBD(design1.df, FALSE)


sim8Tags = simDataChisqPRBD(design = design1.df,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)


sim8Tags$tempS = sim8Tags$tempS[,-2]


edfPlot(sim8Tags, "Within RunBetweenCag:AniResidual")



##########################################################################################
#BIBD Examples
#BIBDdesign4 = optBIBD(nTrt = 4, bRep  = 3, nCag  = 4, tRep  = 2, nPlot = 4, iter = 1000,
#    confoundCag = TRUE)


design.summary.RBD(BIBDdesign4)

(aov.table =summaryAovTwoPhase(BIBDdesign4,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))

sim4Tags = simDataChisqPRBD(design = BIBDdesign4,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)


xx = simDataChisqPRBD(design = BIBDdesign4,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10000,  seed = 1,
            row.MS =rowNames[c(10, 13)])


sim4Tags1$tempS = rbind(sim4Tags$tempS, xx$tempS)
sim4Tags1$tempS = sim4Tags1$tempS[,-2]


edfPlot(sim4Tags1, "Within RunBetweenCag:AniResidual")



#BIBDdesign48 = optBIBD(nTrt = 4, bRep  = 3, nCag  = 4, tRep  = 2, nPlot = 8, iter = 1000,
#    confoundCag = TRUE)

design.summary.RBD(BIBDdesign48)

(aov.table =summaryAovTwoPhase(BIBDdesign48,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")$A)

(rowNames = extractName(rownames(aov.table)))

sim8Tags = simDataChisqPRBD(design = BIBDdesign48,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)
 sfStop()

xx = simDataChisqPRBD(design = BIBDdesign48,
            blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.C = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10000,  seed = 1,
            row.MS =rowNames[c(9, 12)])
  sfStop()

sim8Tags1$tempS = rbind(sim8Tags$tempS, xx$tempS)
sim8Tags1$tempS = sim8Tags1$tempS[,-2]


edfPlot(sim8Tags1, "Within RunBetweenCag:AniResidual")


edfPlotCompare(sim4Tags1, sim8Tags1,  compare = c("4Tag", "8Tag"),
method = "REML", rowNames = "Within RunBetweenCag:AniResidual")
