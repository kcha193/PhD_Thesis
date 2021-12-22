
 
simDataChisqP =
function(design, blk.str2, blk.str1, trt.str, gamma.R,  gamma.A, 
                  VC.resid, nS, nVc, nSim, seed, row.MS = NA, neg.VC = TRUE, checkFun = FALSE){
  
  sfInit(parallel = TRUE, cpus = 4)

  aov.table = summaryAovTwoPhase(design.df = design, blk.str1 =  blk.str1,
                          blk.str2 = blk.str2, trt.str = trt.str, response = rnorm(nrow(design)))
   
     tempS = data.frame()

 
# create progress bar
  for(i in 1:length(gamma.R)){
    gamma.run = gamma.R[i]

    cat("gamma.run = ", gamma.run, "\n")
          
    
   
     test = function(j) {
  
      gamma.ani = gamma.A[j]

      REML.EDF = numeric(nS); LC.EDF = numeric(nS); REAL.EDF  = numeric(nS);
 
      real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.run * VC.resid))  
     index =  1:length(real.VC)
     
      if(any(is.na(real.VC))){
         index = index[-which(is.na(real.VC))]
         real.VC =  real.VC[-which(is.na(real.VC))]
             
      }
 
      set.seed(seed)
      simMS = suppressWarnings(apply(aov.table$A, 1, 
          function(x) sum(fracToNum(x[2:(nVc + 1)]) * real.VC) * rchisq(nSim, fracToNum(x[1]))/fracToNum(x[1])))

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
          
       return(cbind(gamma.run,  gamma.ani, REML.EDF,  LC.EDF, REAL.EDF)["Within RunBetweenAniResidual",]) 
      
    }
     if(checkFun) browser()
 
    sfExport("nSim", "gamma.A", "VC.resid", "gamma.run", "nS", "fracToNum", "nVc", 
    "getVcEDF", "extractName", "row.MS", "neg.VC", "aov.table", "seed")
    
    sfLibrary(MASS)
    
     tempS = rbind(tempS, t(sfSapply(1:length(gamma.A), fun = test)))

  }
   sfStop()

    tempS$gamma.run = as.factor(tempS$gamma.run)
    tempS$gamma.ani  = as.factor(tempS$gamma.ani)
    rownames(tempS) = paste(rep("Within RunBetweenAniResidual", nrow(tempS)), 1:nrow(tempS), sep = "")
    return(list(tempS = tempS))     
  }


design.summary.CRD(design.df, FALSE)

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
 
edfPlot(sim4Tags)



design.summary.CRD(design1.df, FALSE)

old = proc.time()  
sim8Tags = simDataChisqP(design =  design1.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1, neg.VC = TRUE, checkFun = FALSE)
proc.time() - old 
  sfStop()



(aov.table = summaryAovTwoPhase(design1.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(design1.df)))$A)

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

edfPlotCompare(sim4Tags, sim8Tags,  compare = c("4Tag", "8Tag"), method = "REML")

edfPlotCompare(sim4Tags, sim8Tags,  compare = c("4Tag", "8Tag"), method = "LC")
