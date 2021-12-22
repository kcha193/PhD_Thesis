

simData1=
function(design, blk.str2, blk.str1, trt.str,  gamma.A, gamma.P, VC.resid, nS, nVc, nSim, row.MS = NA, neg.VC = TRUE, ...){

  tempS = data.frame()
  tempV = data.frame()

# create progress bar
  for(i in 1:length(gamma.A)){
    gamma.array = gamma.A[i]

    cat("gamma.array = ", gamma.array, "\n")
    for(j in 1:length(gamma.P)){

      gamma.plant = gamma.P[j]

      #REML.EDF = numeric(nS); LC.EDF   = numeric(nS); TRUE.EDF  = numeric(nS);

      REML.VC = numeric(nVc); LC.VC   = numeric(nVc); TRUE.VC  = numeric(nVc); ASREML.VC = numeric(nVc);
      
      y.data = numeric(16)
       cat("gamma.plant = ", gamma.plant, "\n")
       pb <- txtProgressBar(min = 0, max = nSim, style = 3)
       counter = 0
       #simulation
      #while( counter < nSim){
       for( i in 1:nSim){
        setTxtProgressBar(pb, counter)
        SA = with(design, interaction(Set, Array))
        array.eff = rnorm(nlevels(SA), mean = 0, sd = sqrt(gamma.array * VC.resid))
        plant.eff = rnorm(nlevels(design$Plant), mean = 0, sd = sqrt(gamma.plant * VC.resid))
        trt.eff = rep(0, nlevels(design$Treat)) #runif(nlevels(design$Treat), 0, 2)
        dye.eff = rep(0,nlevels(design$Dye)) # runif(nlevels(design$Dye), 0, 1)
        res.eff = rnorm(nrow(design), mean = 0, sd = sqrt(VC.resid))
        gm = 0

        true.VC = c(VC.resid, (gamma.plant * VC.resid),(gamma.array * VC.resid))

        y = gm + with(design, array.eff[SA] + plant.eff[Plant] + dye.eff[Dye] + trt.eff[Treat]) + res.eff

        #XT = matrix(0, ncol = nlevels(design$Treat), nrow = length(y))
        #XT[cbind(1:length(y), design$Treat)] = 1
        
        #XD = matrix(0, ncol = nlevels(design$Dye), nrow = length(y))
        #XD[cbind(1:length(y), design$Dye)] = 1
 
        #X = cbind(1,XT,XD)
        
        #y = (diag(length(y)) - projMat(X)) %*% y
                          
        aov.table = summary.aov.twoPhase(design.df = design, blk.str1 =  blk.str1,
        blk.str2 = blk.str2, trt.str = trt.str, response = y, var.comp = c("Plant", "Set:Array"))

        tmp = try(getVcEDF(aov.table = aov.table, row.MS = row.MS, true.VC = true.VC, neg.VC = neg.VC), TRUE)
        #tmp
        if(class(tmp) =="try-error"){
          #print("ok")
          #k = k - 1
           y.data = rbind(y, y.data) 
          next
        }

       # asreml.fit = asreml(y ~ Treat + Dye, random = ~SA + Plant, data = design, trace = FALSE)
       # ASREML.VC   =  rbind(ASREML.VC,     summary(asreml.fit, all =TRUE)$v[,"component"] )


        
        #tmpS = tmp$S
        tmpV = tmp$V

        #REML.EDF  =  rbind(REML.EDF,     tmpS[,"REML.EDF"] )
        #LC.EDF    =  rbind(LC.EDF  ,     tmpS[,"LC.EDF"]   )
        #TRUE.EDF  =  rbind(TRUE.EDF,     tmpS[,"TRUE.EDF"] )

        REML.VC   =  rbind(REML.VC,      tmpV[,"REML.var.comp"]  )
        LC.VC     =  rbind(LC.VC  ,      tmpV[,"LC.var.comp"]    )
        TRUE.VC   =  rbind(TRUE.VC,      tmpV[,"TRUE.var.comp"]  )

        #y.data   =  rbind(y.data ,       y.data)
        
        counter = counter +1

      }
       close(pb)


       #REML.EDF = apply(REML.EDF[-1,], 2, median)
       #LC.EDF = apply(LC.EDF[-1,], 2, median)
       #TRUE.EDF = apply(TRUE.EDF[-1,], 2, median)
       
       #REML.EDF = REML.EDF[-1,]
       #LC.EDF = LC.EDF[-1,]
       #TRUE.EDF = TRUE.EDF[-1,]
       
       #REML.VC = apply(REML.VC[-1,], 2, mean)
       #LC.VC = apply(LC.VC[-1,], 2, mean)
       #TRUE.VC = apply(TRUE.VC[-1,], 2, mean)

       REML.VC = REML.VC[-1,]
       LC.VC = LC.VC[-1,]
       TRUE.VC = TRUE.VC[-1,]
       
       #tempS = rbind(tempS, cbind(gamma.array, gamma.plant, REML.EDF, LC.EDF , TRUE.EDF) )
       tempV = rbind(tempV, cbind(gamma.array, gamma.plant, REML.VC, LC.VC , TRUE.VC) )
       #fail.y = rbind(tempV, cbind(gamma.array, gamma.plant, y.data) )

    }

  }

    #tempS$gamma.array = as.factor(tempS$gamma.array)
    #tempS$gamma.plant  = as.factor(tempS$gamma.plant)

    tempV$gamma.array = as.factor(tempV$gamma.array)
    tempV$gamma.plant  = as.factor(tempV$gamma.plant)

    #fail.y$gamma.array = as.factor(fail.y$gamma.array)
    #fail.y$gamma.plant  = as.factor(fail.y$gamma.plant)

    #return(list(tempS = tempS, tempV = tempV))
     return(list(tempV = tempV ) )
  }



  plotEDF1 =
function(sim, rowNames = "WithinBetweenPlantResidual", ...){
  require(reshape)

  tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)),]

  mm = melt(tempS.bwAni)

  new.temp = mm[which(!is.na(match(mm$variable, c("REML.EDF", "LC.EDF", "TRUE.EDF")))),]

  new.temp$value = round(new.temp$value, digits = 4)
  colnames(new.temp)[4] = "EDF"


   suppressWarnings(trellis.par.set(layout.heights = list(strip = 1.25),
          superpose.line = list(col=c("black", "red", "blue"), lwd=2, lty=1:3)))

   interval = as.numeric(levels(new.temp$gamma.plant)[seq(1,nlevels(new.temp$gamma.plant),2)])

   new.temp$gamma.plant = as.numeric(as.character(new.temp$gamma.plant))

   xyplot(EDF ~ gamma.plant|gamma.array, group = variable, data = new.temp,lty=1:3, type="l",
   col=c("black", "red", "blue"), lwd=2,
   scales=list(x=list(at=interval, labels=interval, log=TRUE)),
    xlab=expression(paste(sigma[P]^2," / ",sigma^2)),
      key  =
      simpleKey(text = paste(c("Method:  ","", ""), levels(new.temp$variable)),
          points = FALSE, lines = TRUE, cex = 1, columns=3),
     strip = function(factor.levels, par.strip.text, ...)
    {strip.default(factor.levels =
    sapply(levels(new.temp$gamma.array),
    function(x) as.expression(substitute(list(sigma[A]^2/sigma^2==x), list(x=x) ))),
  par.strip.text = trellis.par.get("layout.heights"), ...)},   ... )


}



adjusted =
function(x){
sort(x)[1:ceiling(length(x)*.99)]
}

histVCplots = 
function(sim5, index, r, adjusted = FALSE){
  print(unique(index)[r])
  print("TRUE VC:")
  print(apply(sim5[which(index ==unique(index)[r]),9:11],2,mean))
  
  opar = par(mfrow = c(2,3))
  if(adjusted){
    
    hist(adjusted(sim5[which(index ==unique(index)[r]),3]), 
      main = paste("REML estimates, mean = ", round(mean(adjusted(sim5[which(index ==unique(index)[r]),3])), digits = 4), 
                                    ", median = ", round(median(adjusted(sim5[which(index ==unique(index)[r]),3])), digits = 4), sep = "") , 
      xlab = expression(sigma^2), breaks = 100)
    hist(adjusted(sim5[which(index ==unique(index)[r]),4]), 
      main = paste("REML estimates, mean = ", round(mean(adjusted(sim5[which(index ==unique(index)[r]),4])), digits = 4), 
                                    ", median = ", round(median(adjusted(sim5[which(index ==unique(index)[r]),4])), digits = 4), sep = "") , 
      xlab = expression(sigma[P]^2), breaks = 100)
    hist(adjusted(sim5[which(index ==unique(index)[r]),5]), 
      main = paste("REML estimates, mean = ", round(mean(adjusted(sim5[which(index ==unique(index)[r]),5])), digits = 4), 
                                    ", median = ", round(median(adjusted(sim5[which(index ==unique(index)[r]),5])), digits = 4), sep = "") , 
      xlab = expression(sigma[A]^2), breaks = 100)
    
    hist(adjusted(sim5[which(index ==unique(index)[r]),6]), 
      main = paste("LC estimates, mean = ", round(mean(adjusted(sim5[which(index ==unique(index)[r]),6])), digits = 4), 
                                    ", median = ", round(median(adjusted(sim5[which(index ==unique(index)[r]),6])), digits = 4), sep = "") ,  
      xlab = expression(sigma^2), breaks = 100)
    hist(adjusted(sim5[which(index ==unique(index)[r]),7]), 
      main =  paste("LC estimates, mean = ", round(mean(adjusted(sim5[which(index ==unique(index)[r]),7])), digits = 4), 
                                    ", median = ", round(median(adjusted(sim5[which(index ==unique(index)[r]),7])), digits = 4), sep = "") ,  
      xlab = expression(sigma[P]^2), breaks = 100)
    hist(adjusted(sim5[which(index ==unique(index)[r]),8]), 
      main =  paste("LC estimates, mean = ", round(mean(adjusted(sim5[which(index ==unique(index)[r]),8])), digits = 4), 
                                    ", median = ", round(median(adjusted(sim5[which(index ==unique(index)[r]),8])), digits = 4), sep = "") ,  
      xlab = expression(sigma[A]^2), breaks = 100)

  
  }else{
  
  
    hist(sim5[which(index ==unique(index)[r]),3], 
      main = paste("REML estimates, mean = ", round(mean(sim5[which(index ==unique(index)[r]),3]), digits = 4), 
                                    ", median = ", round(median(sim5[which(index ==unique(index)[r]),3]), digits = 4), sep = "") , 
      xlab = expression(sigma^2), breaks = 100)
    hist(sim5[which(index ==unique(index)[r]),4], 
      main = paste("REML estimates, mean = ", round(mean(sim5[which(index ==unique(index)[r]),4]), digits = 4), 
                                    ", median = ", round(median(sim5[which(index ==unique(index)[r]),4]), digits = 4), sep = "") , 
      xlab = expression(sigma[P]^2), breaks = 100)
    hist(sim5[which(index ==unique(index)[r]),5], 
      main = paste("REML estimates, mean = ", round(mean(sim5[which(index ==unique(index)[r]),5]), digits = 4), 
                                    ", median = ", round(median(sim5[which(index ==unique(index)[r]),5]), digits = 4), sep = "") , 
      xlab = expression(sigma[A]^2), breaks = 100)
    
    hist(sim5[which(index ==unique(index)[r]),6], 
      main = paste("LC estimates, mean = ", round(mean(sim5[which(index ==unique(index)[r]),6]), digits = 4), 
                                    ", median = ", round(median(sim5[which(index ==unique(index)[r]),6]), digits = 4), sep = "") , 
      xlab = expression(sigma^2), breaks = 100)
    hist(sim5[which(index ==unique(index)[r]),7], 
      main =  paste("LC estimates, mean = ", round(mean(sim5[which(index ==unique(index)[r]),7]), digits = 4), 
                                    ", median = ", round(median(sim5[which(index ==unique(index)[r]),7]), digits = 4), sep = "") , 
      xlab = expression(sigma[P]^2), breaks = 100)
    hist(sim5[which(index ==unique(index)[r]),8], 
      main =  paste("LC estimates, mean = ", round(mean(sim5[which(index ==unique(index)[r]),8]), digits = 4), 
                                    ", median = ", round(median(sim5[which(index ==unique(index)[r]),8]), digits = 4), sep = "") ,  
      xlab = expression(sigma[A]^2), breaks = 100)
  
  
  }
  par(opar)
}
