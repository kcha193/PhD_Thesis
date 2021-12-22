



aov.table = summaryAovOnePhase(design.df = design5, blk.str =  "Cag/Ani",
            trt.str = "Trt", response = rnorm(nrow(design5)))

nSim = 10000

gamma.cag = 0.01
gamma.ani = 1000
VC.resid = 1

true.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid))

MS = suppressWarnings(apply(aov.table$A, 1,
          function(x) sum(fracToNum(x[2:4]) * true.VC) * rchisq(nSim, fracToNum(x[1]))/fracToNum(x[1])))


old.time = proc.time()

neg.VC = FALSE

REML.EDF = numeric(3); LC.EDF = numeric(3); REAL.EDF  = numeric(3);
REML.VC = numeric(3); LC.VC = numeric(3); REAL.VC  = numeric(3);

 pb <- txtProgressBar(min = 0, max = nSim, style = 3)
 counter = 0
 #simulation
 while( counter < nSim){
  setTxtProgressBar(pb, counter+1)


  #aov.table$ANOVA[,"MS"] =  MS[counter+1, ]

  tmp = try(getVcEDF(aov.table = aov.table, MS = MS[counter+1, ], true.VC = true.VC, neg.VC = neg.VC), TRUE)

  #print(tmp)
  #tmp
  if(class(tmp) =="try-error"){
  #stop()
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
proc.time() - old.time

summaryEDF = cbind(apply(REAL.EDF[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE)))

colnames(summaryEDF) = c("True.EDF", "REML.EDF.median", "LC.EDF.median",
                                     "REML.EDF.mean", "LC.EDF.mean")
summaryEDF

summaryVC = cbind(apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) median(x, na.rm = TRUE)))

colnames(summaryVC) = c("True.EDF", "REML.VC.mean", "LC.VC.mean",
                                     "REML.VC.median", "LC.VC.median")

summaryVC


#####################################################################################
design = design.df

design$ani = with(design, interaction(Cag, Ani))


REML.EDF = numeric(3); LC.EDF = numeric(3); REAL.EDF  = numeric(3);
REML.VC = numeric(3); LC.VC = numeric(3); REAL.VC  = numeric(3);

old.time = proc.time()
gamma.cag = 0
gamma.ani =  1000
VC.resid = .0001

  real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid))

neg.VC = FALSE


nSim = 10000

cat("gamma.ani = ", gamma.ani, "\n")
pb <- txtProgressBar(min = 0, max = nSim, style = 3)
counter = 0
#simulation
while( counter < nSim){
  setTxtProgressBar(pb, counter+1)
   #set.seed(seed + counter)
  cag.eff = rnorm(nlevels(design$Cag), mean = 0, sd = sqrt(gamma.cag * VC.resid))
  ani.eff = rnorm(nlevels(design$ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
  trt.eff = runif(nlevels(design$Trt), 0, 10)

  res.eff = rnorm(nrow(design), mean = 0, sd = sqrt(VC.resid))
  gm = 10


  y = gm + with(design, cag.eff[Cag] + ani.eff[ani] +  trt.eff[Trt]) + res.eff

  #save.y[counter+1, ] = y

  #y = save.y[counter+1, ]

  aov.table =  summaryAovOnePhase(design.df = design.df, blk.str =  "Cag/Ani",
            trt.str = "Trt", response = y)

  tmp = try(getVcEDF(aov.table = aov.table, row.MS = NA, true.VC = real.VC, neg.VC = neg.VC), TRUE)

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
  proc.time() - old.time

summaryEDF = cbind(apply(REAL.EDF[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE)))

colnames(summaryEDF) = c("True.EDF", "REML.EDF.median", "LC.EDF.median",
                                     "REML.EDF.mean", "LC.EDF.mean")
summaryEDF

summaryVC = cbind(apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) median(x, na.rm = TRUE)))

colnames(summaryVC) = c("True.EDF", "REML.VC.mean", "LC.VC.mean",
                                     "REML.VC.median", "LC.VC.median")
summaryVC

##################################################################################



aov.table = summaryAovOnePhase(design.df = design.df, blk.str =  "Cag/Ani",
            trt.str = "Trt",  response = rnorm(nrow(design.df)))

gamma.cag = 10
gamma.ani =  5
VC.resid = 1
nSim = 10000

real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid))

MS = suppressWarnings(fractions(apply(aov.table$A, 1,
          function(x) sum(fracToNum(x[2:5]) * real.VC) * rchisq(nSim, fracToNum(x[1]))/fracToNum(x[1]))))


old.time = proc.time()

neg.VC = FALSE

REML.EDF = numeric(3); LC.EDF = numeric(3); REAL.EDF  = numeric(3);
REML.VC = numeric(3); LC.VC = numeric(3); REAL.VC  = numeric(3);

 cat("gamma.ani = ", gamma.ani, "\n")
 pb <- txtProgressBar(min = 0, max = nSim, style = 3)
 counter = 0
 #simulation
 while( counter < nSim){
  setTxtProgressBar(pb, counter+1)


  aov.table$ANOVA[,"MS"] =  MS[counter+1, ]

  tmp = try(getVcEDF(aov.table = aov.table, row.MS = NA, true.VC = real.VC, neg.VC = neg.VC), TRUE)

  #print(tmp)
  #tmp
  if(class(tmp) =="try-error"){
  #stop()
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
proc.time() - old.time


summaryEDF = cbind(apply(REAL.EDF[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE)))

colnames(summaryEDF) = c("True.EDF", "REML.EDF.median", "LC.EDF.median",
                                     "REML.EDF.mean", "LC.EDF.mean")
summaryEDF

summaryVC = cbind(apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) median(x, na.rm = TRUE)))

colnames(summaryVC) = c("True.EDF", "REML.VC.mean", "LC.VC.mean",
                                     "REML.VC.median", "LC.VC.median")
summaryVC

######################################################################################################




aov.table = summaryAovTwoPhase(design.df = design6, blk.str1 =  "Cag/Ani",
                          blk.str2 = "Run", trt.str = "Tag + Trt",
                           response = rnorm(nrow(design6)))



Tag.mat = with(design6, as.matrix(table(1:nrow(design6), Tag)))
Trt.mat = with(design6, as.matrix(table(1:nrow(design6), Trt)))
X = cbind(1, Trt.mat, Tag.mat)
   
  run.eff = rnorm(nlevels(design.df$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
    cag.eff = rnorm(nlevels(design.df$Cag), mean = 0, sd = sqrt(gamma.cag * VC.resid))
    ani.eff = rnorm(nlevels(design.df$Ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
    trt.eff = c(1, 6, 3.5, 1, 6, 3.5) #runif(nlevels(design.df$Trt), 0, 10)
    tag.eff = runif(nlevels(design.df$Tag), 0, 1)
    res.eff = rnorm(nrow(design.df), mean = 0, sd = sqrt(VC.resid))
    gm = 10
    
    
    y6 = gm + with(design6, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff
 

aov.table = summaryAovTwoPhase(design.df = design6, blk.str1 =  "Cag/Ani",
                          blk.str2 = "Run", trt.str = "Tag + Trt",
                           response = y6 - as.numeric(X %*%ginv( t(X) %*% X) %*% t(X) %*% y6))


gamma.run  = 0
gamma.cag = 0.01
gamma.ani =  100
VC.resid = 1
nSim = 10000

real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid), (gamma.run * VC.resid))

MS = suppressWarnings(apply(aov.table$A, 1,
          function(x) sum(fracToNum(x[2:5]) * real.VC) * rchisq(nSim, fracToNum(x[1]))/fracToNum(x[1])))


old.time = proc.time()

neg.VC = FALSE

REML.EDF = numeric(5); LC.EDF = numeric(5); REAL.EDF  = numeric(5);
REML.VC = numeric(4); LC.VC = numeric(4); REAL.VC  = numeric(4);

 cat("gamma.ani = ", gamma.ani, "\n")
 pb <- txtProgressBar(min = 0, max = nSim, style = 3)
 counter = 0
 #simulation
 while( counter < nSim){
  setTxtProgressBar(pb, counter+1)


  aov.table$ANOVA[,"MS"] =  MS[counter+1, ]

  tmp = try(getVcEDF(aov.table = aov.table, row.MS = MS[counter+1, ], true.VC = real.VC, 
  neg.VC = neg.VC), TRUE)

  #print(tmp)
  #tmp
  if(class(tmp) =="try-error"){
  #stop()
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
proc.time() - old.time



summaryEDF = cbind(apply(REAL.EDF[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE)))

colnames(summaryEDF) = c("True.EDF", "REML.EDF.median", "LC.EDF.median",
                                     "REML.EDF.mean", "LC.EDF.mean")
summaryEDF

summaryVC = cbind(apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) median(x, na.rm = TRUE)))

colnames(summaryVC) = c("True.EDF", "REML.VC.mean", "LC.VC.mean",
                                     "REML.VC.median", "LC.VC.median")

summaryVC




