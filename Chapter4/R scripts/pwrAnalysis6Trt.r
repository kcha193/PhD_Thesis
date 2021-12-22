
contrCompare1 = function(design1.df, contr){
  n = nrow(design1.df)
  Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))
   
  
  pairMat = contr
  
  
  nBlk = nlevels(design1.df$Run)
  nPlot = nlevels(design1.df$Tag)
  nCag = nlevels(design1.df$Cag)
  
  nAni = nlevels(design1.df$Ani)
  Run.mat = with(design1.df, as.matrix(table(1:n, Run)))
  Tag.mat = with(design1.df, as.matrix(table(1:n, Tag)))
  Cag.mat = with(design1.df, as.matrix(table(1:n, Cag)))
  Ani.mat = with(design1.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
  Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))
  
  C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
  C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
  cage.Rep = n/nCag
  
  C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
  
  Pb = projMat(Run.mat)
  Pb1 = projMat(Tag.mat)
  
  blk.proj = (identityMat(n) - Pb) %*% (identityMat(n) - Pb1)
  
  info.mat = matMulti(blk.proj, Ani.mat, C.cage)
  
    if (any(abs(info.mat) > 1e-07)) {
        PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)),
            C.cage, Ani.mat)
  
        PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
    } else {
        PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani,
            Ani.mat)
  
    }
  
  contr.e = (apply(pairMat, 2, function(c)  t(c) %*% ginv(t(Trt.mat) %*% Trt.mat) %*% c)/apply(pairMat, 2, function(c)  t(c) %*% ginv(t(Trt.mat) %*% PP1 %*% Trt.mat) %*% c))
  
  
  return(contr.e)
}


pairwiseContr = function(design1.df){

  n = nrow(design1.df)
  Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))
  
  states = lapply(1:ncol(Trt.mat), function(x) t(cbind(x, (x+1):ncol(Trt.mat))))
  states = states[-length(states)]
  
  states= t(matrix(c(states, recursive=TRUE), nrow = 2))
  
  pairMat = matrix(0, nrow = ncol(Trt.mat), ncol = 1)
   
  for( i in 1:nrow(states)){
  
    pairMat = cbind(pairMat, pairMat[,1])
    
    pairMat[states[i,1],i+1] = 0.5
    pairMat[states[i,2],i+1] = -0.5
  
  }
  
  
  return(pairMat[,-1])
}

contrCompare = function(design1.df, contr){
  n = nrow(design1.df)
  Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))
  
  Trt.mat = Trt.mat[,sort(colnames(Trt.mat))]
  
  states = lapply(1:ncol(Trt.mat), function(x) t(cbind(x, (x+1):ncol(Trt.mat))))
  states = states[-length(states)]
  
  states= t(matrix(c(states, recursive=TRUE), nrow = 2))
  
   
  pairMat = contr
  
  
  nBlk = nlevels(design1.df$Run)
  nPlot = nlevels(design1.df$Tag)
  nCag = nlevels(design1.df$Cag)
  
  nAni = nlevels(design1.df$Ani)
  Run.mat = with(design1.df, as.matrix(table(1:n, Run)))
  Tag.mat = with(design1.df, as.matrix(table(1:n, Tag)))
  Cag.mat = with(design1.df, as.matrix(table(1:n, Cag)))
  Ani.mat = with(design1.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
  
  C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
  C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
  cage.Rep = n/nCag
  
  C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
  
  Pb = projMat(Run.mat)
  Pb1 = projMat(Tag.mat)
  
  blk.proj = (identityMat(n) - Pb) %*% (identityMat(n) - Pb1)
  
  info.mat = matMulti(blk.proj, Ani.mat, C.cage)
  
    if (any(abs(info.mat) > 1e-07)) {
        PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)),
            C.cage, Ani.mat)
  
        PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
    } else {
        PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani,
            Ani.mat)
  
    }
  
  contr.e = (apply(pairMat, 2, function(c)  t(c) %*% ginv(t(Trt.mat) %*% Trt.mat) %*% c)/apply(pairMat, 2, function(c)  t(c) %*% ginv(t(Trt.mat) %*% PP1 %*% Trt.mat) %*% c))
  
   names(contr.e) =   apply(states, 1, function(x) paste(colnames(Trt.mat)[x], sep = "", collapse = "-"))

  return(contr.e)
}




trtProj = function(design1.df){

n = nrow(design1.df)
Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))



nBlk = nlevels(design1.df$Run)
nPlot = nlevels(design1.df$Tag)
nCag = nlevels(design1.df$Cag)

nAni = nlevels(design1.df$Ani)
Run.mat = with(design1.df, as.matrix(table(1:n, Run)))
Tag.mat = with(design1.df, as.matrix(table(1:n, Tag)))
Cag.mat = with(design1.df, as.matrix(table(1:n, Cag)))
Ani.mat = with(design1.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))

C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
cage.Rep = n/nCag

C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))

Pb = projMat(Run.mat)
Pb1 = projMat(Tag.mat)

blk.proj = (identityMat(n) - Pb) %*% (identityMat(n) - Pb1)

info.mat = matMulti(blk.proj, Ani.mat, C.cage)

  if (any(abs(info.mat) > 1e-07)) {
      PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)),
          C.cage, Ani.mat)

      PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
  } else {
      PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani,
          Ani.mat)

  }

return( ginv(t(Trt.mat) %*% PP1 %*% Trt.mat) %*% t(Trt.mat) %*% PP1)
}


pairContr = data.frame(Design = as.factor(rep(c("6DF", "5DF"), each = 15)), 
                  Eff = c(contrCompare(design6, pairwiseContr(design6)),  
                  contrCompare(design5, pairwiseContr(design5))))

boxplot(Eff ~ Design, pairContr)

qplot( Design,Eff, data = pairContr,
    main = "Boxplot of the efficiency factors of all pairwise contrasts", 
    xlab =  "Design", ylab = "Efficiency factor", geom = "boxplot")
                  
p  + geom_boxplot( )

###################################################################################

(anovaTable6 = summaryAovTwoPhase(design.df = design6 ,blk.str2 = "Run",
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)

(anovaTable5 = summaryAovTwoPhase(design.df = design5 ,blk.str2 = "Run",
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)

temp6 = data.frame()
temp5 = data.frame()

nSim = 10000

VC.resid = 1

gamma.run =  10
gamma.cag =  5

gamma.A = c(10^((-4:4)/2))
trt.Eff = rbind(c(2,2,0,0,2,0),
                c(6, 1, 6, 6, 11, 6),
                c(3, 1, 3, 5, 3, 3),
                c(6, 1, 6, 11, 6, 6),
                 c(1, 1, 1, 1, 1, 1),
                c(1, 2, 3, 4, 5, 6),
                c(1, 5, 9, 13, 17, 21))

for(i in 1:1){
cat("trt.Eff ", i, "\n")

for(j in 1:length(gamma.A)){

gamma.ani =  gamma.A[j]


# create progress bar
MS6 = rep(0, 2)
MS5 = rep(0, 2)

save.y6 = rep(0, length(y6))
save.y5 = rep(0, length(y6))


cat("gamma.ani = ", gamma.ani, "\n")
pb <- txtProgressBar(min = 0, max = nSim, style = 3)
counter = 0
# simulation

while (counter < nSim) {
    setTxtProgressBar(pb, counter + 1)

    run.eff = rnorm(nlevels(design.df$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
    cag.eff = rnorm(nlevels(design.df$Cag), mean = 0, sd = sqrt(gamma.cag * VC.resid))
    ani.eff = rnorm(nlevels(design.df$Ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
    trt.eff = trt.Eff[i,] #runif(nlevels(design.df$Trt), 0, 10)
    tag.eff = runif(nlevels(design.df$Tag), 0, 1)
    res.eff = rnorm(nrow(design.df), mean = 0, sd = sqrt(VC.resid))
    gm = 10


    y6 = gm + with(design6, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] +
            tag.eff[Tag] + trt.eff[Trt]) + res.eff

    save.y6 =  cbind(save.y6, y6)

    y5 = gm + with(design5, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] +
            tag.eff[Tag] + trt.eff[Trt]) + res.eff

    save.y5 =  cbind(save.y5, y5)


    aov.table6 = summaryAovTwoPhase(design.df = design6, blk.str2 = "Run",
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt",
      response = y6)

    aov.table5 = summaryAovTwoPhase(design.df = design5, blk.str2 = "Run",
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt",
      response = y5)




    MS6 = cbind(MS6, as.numeric(aov.table6$A[10:11, "MS"]))
    MS5 = cbind(MS5, as.numeric(aov.table5$A[12:13, "MS"]))

    counter = counter + 1
}
close(pb)


   MS6 = MS6[, -1]
        MS5 = MS5[, -1]
        
        save.y6 = save.y6[, -1]
        save.y5 = save.y5[, -1]
        
        
        trtVar = sapply(MS6[2, ], function(x) sqrt(2 * x/(6 * contr6.e)))
        
        
        trt.Est6 = apply(save.y6, 2, function(x) 2 * as.numeric(trtEff6 %*% x) %*% pairwiseContr(design6))
               
        
        p.value6 = pt(-abs(trt.Est6/trtVar), df = 6) * 2
        
        power6 = apply(p.value6, 1, function(x) sum(x < 0.05)/length(x))
        
        
        trtVar = sapply(MS5[2, ], function(x) sqrt(2 * x/(6 * contr5.e)))
        
        
        trt.Est5 = apply(save.y5, 2, function(x) 2 * as.numeric(trtEff5 %*% x) %*% pairwiseContr(design5))
             
        p.value5 = pt(-abs(trt.Est5/trtVar), df = 5) * 2
        
        power5 = apply(p.value5, 1, function(x) sum(x < 0.05)/length(x))
        
        #data.frame(contr = apply(states, 1, function(x) paste(letters[x],  sep = "", collapse = "-")), power6, power5)
        

temp6 = rbind(temp6, c(i,
              gamma.ani,
              sum((1 - pf(MS6[1, ]/MS6[2, ], 5, 6)) < 0.05)/ nSim, power6))

temp5 = rbind(temp5, c(i,
              gamma.ani, sum((1 - pf(MS5[1, ]/MS5[2, ], 5, 5)) < 0.05)/nSim, power5))

        }

    }

colnames(temp6) = c("Trt.Eff", "gamma.ani", "F-test power",
  apply(states, 1, function(x) paste(x, sep = "", collapse = "-")))

colnames(temp5) = c("Trt.Eff", "gamma.ani", "F-test power",
  apply(states, 1, function(x) paste(x, sep = "", collapse = "-")))

 eff = round(c(0, 0, 1/mean(1/contrCompare(design6, pairwiseContr(design6))), contrCompare(design6, pairwiseContr(design6))), 3)

temp6 = rbind(eff, temp6)

  
  eff =  round(c(0, 0, 1/mean(1/contrCompare(design5, pairwiseContr(design5))), contrCompare(design5, pairwiseContr(design5))), 3)

temp5 = rbind(eff, temp5)

 
temp6[,"Trt.Eff"] = c(0,apply(trt.Eff[temp6[,"Trt.Eff"][-1],], 1, function(x) paste(x, sep = "", collapse = ",")))

temp5[,"Trt.Eff"] = c(0,apply(trt.Eff[temp5[,"Trt.Eff"][-1],], 1, function(x) paste(x, sep = "", collapse = ",")))



##########################################################################################
#Reference

ref.temp6 = data.frame()
ref.temp5 = data.frame()

    sig.level = 0.05

gamma.A = c(10^((-4:4)/2))
trt.Eff = rbind(c(1, 5, 3, 3, 3, 3),
                c(1, 11, 6, 6, 6, 6),
                c(1, 1, 1, 5, 5, 5),
                c(1, 1, 1, 11, 11, 11),
                c(1, 2, 3, 4, 5, 6),
                c(1, 5, 9, 13, 17, 21))

 for(i in 1:nrow(trt.Eff)){
cat("trt.Eff ", i, "\n")
 trt.eff =  trt.Eff[i,]

for(j in 1:length(gamma.A)){

gamma.ani =  gamma.A[j]

   groups = length(trt.eff)
    between.var = var( trt.eff)
    within.var = 1 + 2*gamma.ani

    lambda <- (groups  -1 ) *  6 * 2600/3169 * (between.var/within.var)
      lower.df = 6

     #lambda <- (groups - 1) *  MS6[10,] /(MS6[11, ]/(2600/3169))

fPower6 =     pf(qf(sig.level, groups - 1, lower.df , lower.tail = FALSE),
            groups - 1, lower.df , lambda, lower.tail = FALSE)

 contr.e = contrCompare(design6, pairwiseContr(design6))

trtVar = sqrt(2*(1 + 2*gamma.ani)/(6* contr.e))


 trt.Est6 =   2 * as.numeric(trtProj(design6) %*% trt.eff[design6$Trt]) %*% pairwiseContr(design6)

tPower6 = pt(qt(sig.level/2, lower.df, lower.tail = FALSE), lower.df, 
ncp =  abs(trt.Est6) /trtVar, lower.tail = FALSE)


      lower.df = 5
      lambda <- (groups - 1)  * 6  *  13160/15859 * (between.var/within.var)

fPower5 =  pf(qf(sig.level, groups - 1, lower.df , lower.tail = FALSE),
            groups - 1, lower.df , lambda, lower.tail = FALSE)


   contr.e = contrCompare(design5, pairwiseContr(design5))

trtVar = sqrt(2*(1 + 2*gamma.ani)/(6* contr.e))


 trt.Est5 =  2 * as.numeric(trtProj(design5) %*% trt.eff[design5$Trt]) %*% pairwiseContr(design6)

tPower5 =pt(qt(sig.level/2, lower.df, lower.tail = FALSE), lower.df, ncp = abs(trt.Est5) /trtVar, lower.tail = FALSE)

ref.temp6 = rbind(ref.temp6, round(c(i, gamma.ani, fPower6, tPower6), 3))

ref.temp5 = rbind(ref.temp5,  round(c(i, gamma.ani, fPower5, tPower5), 3))
}

}            

colnames(ref.temp6) = c("Trt.Eff", "gamma.ani", "F-test power",
  apply(states, 1, function(x) paste(x, sep = "", collapse = "-")))

colnames(ref.temp5) = c("Trt.Eff", "gamma.ani", "F-test power",
  apply(states, 1, function(x) paste(x, sep = "", collapse = "-")))

 eff = round(c(0, 0, 1/mean(1/contrCompare1(design6, basicContr(design6))), contrCompare(design6, pairwiseContr(design6))), 4)

ref.temp6 = rbind(eff, ref.temp6)


  eff =  round(c(0, 0, 1/mean(1/contrCompare1(design5, basicContr(design5))), contrCompare(design5, pairwiseContr(design5))), 4)

ref.temp5 = rbind(eff, ref.temp5)


ref.temp6[,"Trt.Eff"] = c(0,apply(trt.Eff[ref.temp6[,"Trt.Eff"][-1],], 1, function(x) paste(x, sep = "", collapse = ",")))

ref.temp5[,"Trt.Eff"] = c(0,apply(trt.Eff[ref.temp5[,"Trt.Eff"][-1],], 1, function(x) paste(x, sep = "", collapse = ",")))

ref.temp6
ref.temp5

##################################################################################
library(snowfall)   #Use parallel programming to speed up the computation by 4 times
 sfInit( parallel=TRUE, cpus=4 )

(anovaTable6 = summaryAovTwoPhase(design.df = design6 ,blk.str2 = "Run",
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)

(anovaTable5 = summaryAovTwoPhase(design.df = design5 ,blk.str2 = "Run",
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)

levels(design5$Trt) = c("a", "f", "b", "d", "e", "c")

  contr6.e = contrCompare(design6, pairwiseContr(design6))
   contr5.e = contrCompare(design5, pairwiseContr(design5))

     trtEff6 = trtProj(design6)
   trtEff5 = trtProj(design5)

temp6 = data.frame()
temp5 = data.frame()

nSim = 10000

VC.resid = 1

gamma.run =  10
gamma.cag =  5

gamma.A = c(10^((-4:4)/2))
trt.Eff = rbind(c(2,2,2,0,0,0),
                c(6, 1, 6, 6, 11, 6),
                c(3, 1, 3, 5, 3, 3),
                c(6, 1, 6, 11, 6, 6),
                c(1, 2, 3, 4, 5, 6),
                c(1, 5, 9, 13, 17, 21))

for(i in 1:1){
  cat("trt.Eff ", i, "\n")

  trt.eff = trt.Eff[i,]

test = function(j) {
    sourceDir(path = "C:/Users/Kevin/Dropbox/Packaging tools/packaging/infoDecompuTE/R",
    trace = FALSE)


    # create progress bar
    MS6 = rep(0, 2)
    MS5 = rep(0,2)

    save.y6 = rep(0, length(y6))
    save.y5 = rep(0, length(y6))

counter = 0
# simulation

while (counter < nSim) {

        run.eff = rnorm(nlevels(design.df$Run), mean = 0, sd = sqrt(10 * VC.resid))
        cag.eff = rnorm(nlevels(design.df$Cag), mean = 0, sd = sqrt(5 * VC.resid))
        ani.eff = rnorm(nlevels(design.df$Ani), mean = 0, sd = sqrt(gamma.A[j] * VC.resid))
        #trt.eff = trt.Eff[i, ]  #runif(nlevels(design.df$Trt), 0, 10)
        tag.eff = runif(nlevels(design.df$Tag), 0, 1)
        res.eff = rnorm(nrow(design.df), mean = 0, sd = sqrt(VC.resid))
        gm = 10

        y6 = gm + with(design6, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + tag.eff[Tag] +
            trt.eff[Trt]) + res.eff

        save.y6 = cbind(save.y6, y6)

        y5 = gm + with(design5, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + tag.eff[Tag] +
            trt.eff[Trt]) + res.eff

        save.y5 = cbind(save.y5, y5)


        aov.table6 = summaryAovTwoPhase(design.df = design6, blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt", response = y6)

        aov.table5 = summaryAovTwoPhase(design.df = design5, blk.str2 = "Run", blk.str1 = "Cag/Ani",
            trt.str = "Tag + Trt", response = y5)


        MS6 = cbind(MS6, as.numeric(aov.table6$A[10:11, "MS"]))
        MS5 = cbind(MS5, as.numeric(aov.table5$A[12:13, "MS"]))

        counter = counter + 1
    }

   MS6 = MS6[, -1]
        MS5 = MS5[, -1]
        
        save.y6 = save.y6[, -1]
        save.y5 = save.y5[, -1]
        
        
        trtVar = sapply(MS6[2, ], function(x) sqrt(2 * x/(6 * contr6.e)))
        
        
        trt.Est6 = apply(save.y6, 2, function(x) 2 * as.numeric(trtEff6 %*% x) %*% pairwiseContr(design6))
               
        
        p.value6 = pt(-abs(trt.Est6/trtVar), df = 6) * 2
        
        power6 = apply(p.value6, 1, function(x) sum(x < 0.05)/length(x))
        
        
        trtVar = sapply(MS5[2, ], function(x) sqrt(2 * x/(6 * contr5.e)))
        
        
        trt.Est5 = apply(save.y5, 2, function(x) 2 * as.numeric(trtEff5 %*% x) %*% pairwiseContr(design5))
             
        p.value5 = pt(-abs(trt.Est5/trtVar), df = 5) * 2
        
        power5 = apply(p.value5, 1, function(x) sum(x < 0.05)/length(x))
        
        #data.frame(contr = apply(states, 1, function(x) paste(letters[x],  sep = "", collapse = "-")), power6, power5)
        
        return(c(gamma.A[j], sum((1 - pf(MS6[1, ]/MS6[2, ], 5, 6)) < 0.05)/nSim, 
            power6, sum((1 - pf(MS5[1, ]/MS5[2, ], 5, 5)) < 0.05)/nSim, 
            power5))
}

  sfExport("gamma.A", "trt.eff", "VC.resid", "design5", "design6", "sourceDir", "anovaTable6", 
        "anovaTable5", "y6", "nSim", "pairwiseContr", "contr5.e", "contr6.e", "trtProj", "trtEff6", 
        "trtEff5", "design.df")

    temp = t(sfSapply(1:length(gamma.A),  fun = test))

   temp6 = rbind(temp6, cbind(i, temp[, 1:17]))
    
    temp5 = rbind(temp5, cbind(i, temp[, c(1, 18:33)]))

}
    
 sfStop()
 
 states = states[,1:2]


colnames(temp6) = c("Trt.Eff", "gamma.ani", "F-test power", apply(states, 1, function(x) paste(letters[x], 
    sep = "", collapse = "-")))

colnames(temp5) = c("Trt.Eff", "gamma.ani", "F-test power", apply(states, 1, function(x) paste(letters[x], 
    sep = "", collapse = "-")))

eff = round(c(0, 0, 1/mean(1/contrCompare1(design6, basicContr(design6)[,-6])), contrCompare(design6, 
    pairwiseContr(design6))), 5)

temp6 = rbind(eff, temp6)


eff = round(c(0, 0, 1/mean(1/contrCompare1(design5, basicContr(design5)[,-6])), contrCompare(design5, 
    pairwiseContr(design5))), 5)

temp5 = rbind(eff, temp5)


temp6[, "Trt.Eff"] = c(0, apply(trt.Eff[temp6[, "Trt.Eff"][-1], ], 1, function(x) paste(x, sep = "", 
    collapse = ",")))

temp5[, "Trt.Eff"] = c(0, apply(trt.Eff[temp5[, "Trt.Eff"][-1], ], 1, function(x) paste(x, sep = "", 
    collapse = ",")))

write.csv(temp6, "temp65.csv")
write.csv(temp5, "temp55.csv")

qplot(1:30, qt(0.05, df = 1:30, lower.tail = FALSE), geom = "path", 
xlab = "Residual DF", ylab = "Critical T-statistics", 
main = "Crtical T-statistics to yield significant test at 0.05 verus reidual DF")










cbind(contrCompare(design6, pairwiseContr(design6)),
contrCompare(design5, pairwiseContr(design5)))

 levels(newDesign5$Trt) = c("a", "f", "c", "d", "e", "b")
 
contrCompare(newDesign5, pairwiseContr(newDesign5))

