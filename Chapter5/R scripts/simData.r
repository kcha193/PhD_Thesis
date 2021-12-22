
  library(lattice)
  library(reshape)
   library(MASS)

 
simDataChisq =
function(design, blk.str2, blk.str1, trt.str, gamma.R,  gamma.A, 
                  VC.resid, nS, nVc, nSim, seed, row.MS = NA, neg.VC = TRUE, checkFun = FALSE, average = TRUE){

  aov.table = summaryAovTwoPhase(design.df = design, blk.str1 =  blk.str1,
                          blk.str2 = blk.str2, trt.str = trt.str, response = rnorm(nrow(design)))
 
  tempS = data.frame()
  tempV = data.frame()

# create progress bar
  for(i in 1:length(gamma.R)){
    gamma.run = gamma.R[i]

    cat("gamma.run = ", gamma.run, "\n")

    for(j in 1:length(gamma.A)){

      gamma.ani = gamma.A[j]

      REML.EDF = numeric(nS); LC.EDF = numeric(nS); REAL.EDF  = numeric(nS);

      REML.VC = numeric(nVc); LC.VC = numeric(nVc); REAL.VC  = numeric(nVc);

      real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.run * VC.resid))  
      
      if(any(is.na(real.VC))){
         real.VC =  real.VC[-which(is.na(real.VC))]     
      }

  
      set.seed(seed)
      simMS = suppressWarnings(apply(aov.table$A, 1,
                                     function(x)
                                       sum(fracToNum(x[2:(nVc + 1)]) * real.VC) * 
                                       rchisq(nSim, fracToNum(x[1])) /
                                       fracToNum(x[1])))
      
      if(checkFun) browser()

       cat("gamma.ani = ", gamma.ani, "\n")
       pb <- txtProgressBar(min = 0, max = nSim, style = 3)
       counter = 0
       #simulation
       while( counter < nSim){
        setTxtProgressBar(pb, counter+1)      
       
        tmp = try(getVcEDF(aov.table = aov.table, MS = simMS[counter+1, ], true.VC = real.VC, 
                    neg.VC = neg.VC, row.MS = row.MS), TRUE)
        
        print(tmp)
        #tmp
        if(class(tmp) =="try-error"){
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
            if(checkFun) browser()

      if(average){
        REML.EDF = apply(REML.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
       LC.EDF = apply(LC.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
       REAL.EDF = apply(REAL.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
      } else {
      REML.EDF = apply(REML.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
       LC.EDF = apply(LC.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
       REAL.EDF = apply(REAL.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
      }
     
      REML.VC = apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
       LC.VC = apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
       REAL.VC = apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
       
       tempS = rbind(tempS, cbind(gamma.run,  gamma.ani, REML.EDF,  LC.EDF, REAL.EDF) )
       tempV = rbind(tempV, cbind(gamma.run,  gamma.ani,  REML.VC,  LC.VC,  REAL.VC) )

    }
    
  }
  
    tempS$gamma.run = as.factor(tempS$gamma.run)
    tempS$gamma.ani  = as.factor(tempS$gamma.ani)

    tempV$gamma.run = as.factor(tempV$gamma.run)

    tempV$gamma.ani  = as.factor(tempV$gamma.ani)

    return(list(tempS = tempS, tempV = tempV))
  
  }




simDataChisqNew =
  function(design, blk.str2, blk.str1, trt.str, gamma.R,  gamma.A, 
           VC.resid, nS, nVc, nSim, seed, row.MS = NA, neg.VC = TRUE, 
           checkFun = FALSE, average = TRUE, getVcEDF){
    
    aov.table = summaryAovTwoPhase(design.df = design, blk.str1 =  blk.str1,
                                   blk.str2 = blk.str2, trt.str = trt.str, response = rnorm(nrow(design)))
    
    tempS = data.frame()
    tempV = data.frame()
    
    # create progress bar
    for(i in 1:length(gamma.R)){
      gamma.run = gamma.R[i]
      
      cat("gamma.run = ", gamma.run, "\n")
      
      for(j in 1:length(gamma.A)){
        
        gamma.ani = gamma.A[j]
        
        REML.EDF = numeric(nS); LC.EDF = numeric(nS); REAL.EDF  = numeric(nS);
        
        REML.VC = numeric(nVc); LC.VC = numeric(nVc); REAL.VC  = numeric(nVc);
        
        real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.run * VC.resid))  
        
        if(any(is.na(real.VC))){
          real.VC =  real.VC[-which(is.na(real.VC))]     
        }
        
        
        set.seed(seed)
        simMS = suppressWarnings(apply(aov.table$A, 1,
                                       function(x)
                                         sum(fracToNum(x[2:(nVc + 1)]) * real.VC) * 
                                         rchisq(nSim, fracToNum(x[1])) /
                                         fracToNum(x[1])))
        
        if(checkFun) browser()
        
        cat("gamma.ani = ", gamma.ani, "\n")
        pb <- txtProgressBar(min = 0, max = nSim, style = 3)
        counter = 0
        #simulation
        while( counter < nSim){
          setTxtProgressBar(pb, counter+1)      
          
          tmp = try(getVcEDF(aov.table = aov.table, MS = simMS[counter+1, ], true.VC = real.VC, 
                             neg.VC = neg.VC, row.MS = row.MS), TRUE)
          
          print(tmp)
          #tmp
          if(class(tmp) =="try-error"){
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
        if(checkFun) browser()
        
        if(average){
          REML.EDF = apply(REML.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
          LC.EDF = apply(LC.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
          REAL.EDF = apply(REAL.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
        } else {
          REML.EDF = apply(REML.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
          LC.EDF = apply(LC.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
          REAL.EDF = apply(REAL.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
        }
        
        REML.VC = apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
        LC.VC = apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
        REAL.VC = apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
        
        tempS = rbind(tempS, cbind(gamma.run,  gamma.ani, REML.EDF,  LC.EDF, REAL.EDF) )
        tempV = rbind(tempV, cbind(gamma.run,  gamma.ani,  REML.VC,  LC.VC,  REAL.VC) )
        
      }
      
    }
    
    tempS$gamma.run = as.factor(tempS$gamma.run)
    tempS$gamma.ani  = as.factor(tempS$gamma.ani)
    
    tempV$gamma.run = as.factor(tempV$gamma.run)
    
    tempV$gamma.ani  = as.factor(tempV$gamma.ani)
    
    return(list(tempS = tempS, tempV = tempV))
    
  }

#######################################################################################

simDataChisq1 =
function(design, blk.str2, blk.str1, trt.str, gamma.R,  gamma.C, gamma.A, 
                  VC.resid, nS, nVc, nSim, seed, row.MS = NA, neg.VC = TRUE, checkFun = FALSE){

  aov.table = summaryAovTwoPhase(design.df = design, blk.str1 =  blk.str1,
                          blk.str2 = blk.str2, trt.str = trt.str, response = rnorm(nrow(design)))
   
  tempS = data.frame()
  tempV = data.frame()

# create progress bar
  for(i in 1:length(gamma.R)){
    gamma.run = gamma.R[i]

    cat("gamma.run = ", gamma.run, "\n")
   for(k in 1:length(gamma.C)){
    gamma.cag = gamma.C[k]
  
     cat("gamma.cag = ", gamma.cag, "\n") 
    for(j in 1:length(gamma.A)){

      gamma.ani = gamma.A[j]

      REML.EDF = numeric(nS); LC.EDF = numeric(nS); REAL.EDF  = numeric(nS);

      REML.VC = numeric(nVc); LC.VC = numeric(nVc); REAL.VC  = numeric(nVc);

      real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid), (gamma.run * VC.resid))  
      index =  1:length(real.VC)
      if(any(is.na(real.VC))){
         index = index[-which(is.na(real.VC))]
         real.VC =  real.VC[-which(is.na(real.VC))]
             
      }

      set.seed(seed)
      simMS = suppressWarnings(apply(aov.table$A, 1, 
          function(x) sum(fracToNum(x[index+1]) * real.VC) * rchisq(nSim, fracToNum(x[1]))/fracToNum(x[1])))
     
     if(checkFun) browser()
       
       cat("gamma.ani = ", gamma.ani, "\n")
       pb <- txtProgressBar(min = 0, max = nSim, style = 3)
       counter = 0
       #simulation
       while( counter < nSim){
        setTxtProgressBar(pb, counter+1)      
       
        tmp = try(getVcEDF(aov.table = aov.table, MS = simMS[counter+1, ], true.VC = real.VC, 
                    neg.VC = neg.VC, row.MS = row.MS), TRUE)
        
        #print(tmp)
        #tmp
        if(class(tmp) =="try-error"){
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
       
        if(checkFun) browser()

      #REML.EDF = apply(REML.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
     #  LC.EDF = apply(LC.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
     #  REAL.EDF = apply(REAL.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
      
           REML.EDF = apply(REML.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
       LC.EDF = apply(LC.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
       REAL.EDF = apply(REAL.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
   
      REML.VC = apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
       LC.VC = apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
       REAL.VC = apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
       
       tempS = rbind(tempS, cbind(gamma.run,  gamma.cag, gamma.ani, REML.EDF,  LC.EDF, REAL.EDF) )
       tempV = rbind(tempV, cbind(gamma.run, gamma.cag, gamma.ani,  REML.VC,  LC.VC,  REAL.VC) )

    }
    }
    
  }
  
    tempS$gamma.run = as.factor(tempS$gamma.run)
      tempS$gamma.cag  = as.factor(tempS$gamma.cag)

    tempS$gamma.ani  = as.factor(tempS$gamma.ani)

    tempV$gamma.run = as.factor(tempV$gamma.run)

      tempS$gamma.cag  = as.factor(tempS$gamma.cag)

    tempV$gamma.ani  = as.factor(tempV$gamma.ani)

    return(list(tempS = tempS, tempV = tempV))
  
  }
 
simData =
function(design, blk.str2, blk.str1, trt.str, gamma.R,  gamma.A, 
                  VC.resid, nS, nVc, nSim, seed, row.MS = NA, neg.VC = TRUE, ...){

  tempS = data.frame()
  tempV = data.frame()

# create progress bar
  for(i in 1:length(gamma.R)){
    gamma.run = gamma.R[i]

    cat("gamma.run = ", gamma.run, "\n")

    for(j in 1:length(gamma.A)){

      gamma.ani = gamma.A[j]

      REML.EDF = numeric(nS); LC.EDF = numeric(nS); REAL.EDF  = numeric(nS);

      REML.VC = numeric(nVc); LC.VC = numeric(nVc); REAL.VC  = numeric(nVc);

       cat("gamma.ani = ", gamma.ani, "\n")
       pb <- txtProgressBar(min = 0, max = nSim, style = 3)
       counter = 0
       #simulation
       while( counter < nSim){
        setTxtProgressBar(pb, counter+1)
         set.seed(seed + counter)
        run.eff = rnorm(nlevels(design$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
        ani.eff = rnorm(nlevels(design$Ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
        trt.eff = runif(nlevels(design$Trt), 0, 2)
        tag.eff = runif(nlevels(design$Tag), 0, 1)
        res.eff = rnorm(nrow(design), mean = 0, sd = sqrt(VC.resid))
        gm = 10

        real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.run * VC.resid))

        y = gm + with(design, run.eff[Run] + ani.eff[Ani] +  tag.eff[Tag] + trt.eff[Trt]) + res.eff

        aov.table = summaryAovTwoPhase(design.df = design, blk.str1 =  blk.str1,
                          blk.str2 = blk.str2, trt.str = trt.str, response = y)

        tmp = try(getVcEDF(aov.table = aov.table, row.MS = NA, true.VC = real.VC, neg.VC = neg.VC), TRUE)
        
        #print(tmp)
        #tmp
        if(class(tmp) =="try-error"){
          counter = counter - 1
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
       
     REML.EDF = apply(REML.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
       LC.EDF = apply(LC.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
       REAL.EDF = apply(REAL.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
      
      #     REML.EDF = apply(REML.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
      # LC.EDF = apply(LC.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
      # REAL.EDF = apply(REAL.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
      
  
  REML.VC = apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
       LC.VC = apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
       REAL.VC = apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
       
       tempS = rbind(tempS, cbind(gamma.run,  gamma.ani, REML.EDF,  LC.EDF, REAL.EDF) )
       tempV = rbind(tempV, cbind(gamma.run,  gamma.ani,  REML.VC,  LC.VC,  REAL.VC) )

    }
    
  }
  
    tempS$gamma.run = as.factor(tempS$gamma.run)
    tempS$gamma.ani  = as.factor(tempS$gamma.ani)

    tempV$gamma.run = as.factor(tempV$gamma.run)

    tempV$gamma.ani  = as.factor(tempV$gamma.ani)

    return(list(tempS = tempS, tempV = tempV))
  
  }


simData1=
function(design, blk.str2, blk.str1, trt.str, gamma.R, gamma.C, gamma.A, 
                  VC.resid, nS, nVc, nSim, row.MS = NA, neg.VC = TRUE, ...){

  tempS = data.frame()
  tempV = data.frame()

# create progress bar
  for(i in 1:length(gamma.R)){
    gamma.run = gamma.R[i]

    cat("gamma.run = ", gamma.run, "\n")
    for(k in 1:length(gamma.C)){
     gamma.cag = gamma.C[k]
     cat("gamma.cag = ", gamma.cag, "\n")

    for(j in 1:length(gamma.A)){

      gamma.ani = gamma.A[j]

      REML.EDF = numeric(nS); LC.EDF = numeric(nS); REAL.EDF  = numeric(nS);

      REML.VC = numeric(nVc); LC.VC = numeric(nVc); REAL.VC  = numeric(nVc);

       cat("gamma.ani = ", gamma.ani, "\n")
       pb <- txtProgressBar(min = 0, max = nSim, style = 3)
       counter = 0
       #simulation
       while( counter < nSim){
        setTxtProgressBar(pb, counter+1)
        design$ani =  with(design, interaction(Cag, Ani))

        run.eff = rnorm(nlevels(design$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
        cag.eff = rnorm(nlevels(design$Cag), mean = 0, sd = sqrt(gamma.cag * VC.resid))
        ani.eff = rnorm(nlevels(design$ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
        trt.eff = runif(nlevels(design$Trt), 0, 2)
        tag.eff = runif(nlevels(design$Tag), 0, 1)
        res.eff = rnorm(nrow(design), mean = 0, sd = sqrt(VC.resid))
        gm = 10

        real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid), (gamma.run * VC.resid))

        y = gm + with(design, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] +  tag.eff[Tag] + trt.eff[Trt]) + res.eff

        aov.table = summaryAovTwoPhase(design.df = design, blk.str1 =  blk.str1,
                          blk.str2 = blk.str2, trt.str = trt.str, response = y)

        tmp = try(getVcEDF(aov.table = aov.table, row.MS = NA, true.VC = real.VC, neg.VC = neg.VC), TRUE)
        #tmp = try(getVcEDF.modified(aov.table = aov.table, row.MS = NA, true.VC = real.VC, neg.VC = neg.VC), TRUE)
        
        #tmp
        
        if(class(tmp) =="try-error"){
          counter = counter - 1
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
       
     REML.EDF = apply(REML.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
       LC.EDF = apply(LC.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
       REAL.EDF = apply(REAL.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
       
       REML.VC = apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
       LC.VC = apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
       REAL.VC = apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))

       tempS = rbind(tempS, cbind(gamma.run, gamma.cag, gamma.ani, REML.EDF, LC.EDF, REAL.EDF) )
       tempV = rbind(tempV, cbind(gamma.run, gamma.cag, gamma.ani, REML.VC, LC.VC, REAL.VC) )

    }
    }
  }
  
    tempS$gamma.run = as.factor(tempS$gamma.run)
    tempS$gamma.cag = as.factor(tempS$gamma.cag)
    tempS$gamma.ani  = as.factor(tempS$gamma.ani)

    tempV$gamma.run = as.factor(tempV$gamma.run)
    tempV$gamma.cag = as.factor(tempV$gamma.cag)
    tempV$gamma.ani  = as.factor(tempV$gamma.ani)

    return(list(tempS = tempS, tempV = tempV))
  
  }



