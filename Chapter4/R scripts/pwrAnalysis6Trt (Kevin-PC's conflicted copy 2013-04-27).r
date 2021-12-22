
(anovaTable6 = summaryAovTwoPhase(design.df = design6 ,blk.str2 = "Run",
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)

(anovaTable5 = summaryAovTwoPhase(design.df = design5 ,blk.str2 = "Run",
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)

temp6 = data.frame()
temp5 = data.frame()

nSim = 1000

VC.resid = 1

gamma.run =  10
gamma.cag =  5

gamma.A = c(10^((-4:4)/2))
#                 a  b  c  d  e  d
trt.Eff = rbind(c(3, 1, 3, 3, 5, 3),
                c(6, 1, 6, 6, 11, 6),
                c(3, 1, 3, 5, 3, 3),
                c(6, 1, 6, 11, 6, 6),
                c(1, 1, 1, 1, 1, 1),
                c(1, 2, 3, 4, 5, 6),
                c(1, 5, 9, 13, 17, 21))

for(i in 1:nrow(trt.Eff)){
cat("trt.Eff ", i, "\n")

for(j in 1:length(gamma.A)){

gamma.ani =  gamma.A[j]


# create progress bar
MS6 = rep(0, nrow(anovaTable6))
MS5 = rep(0, nrow(anovaTable5))

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

save.y6 =   save.y6[, -1]
save.y5 =   save.y5[, -1]


contr.e = contrCompare(design6, pairwiseContr(design6))


trtVar = rowMeans(sapply(MS6[2,], function(x) sqrt(2*x/(6* contr.e))))


 trt.Est6 = apply(save.y6, 2, function(x) 2*tapply(x, design6$Trt, mean) %*% pairwiseContr(design6)* contr.e)

 p.value6 =  apply(trt.Est6, 2, function(x)  pt(-abs(x/trtVar), df = 6) * 2)

 power6 = apply(p.value6, 1, function(x) sum(x<0.05)/length(x))


contr.e = contrCompare(design5, pairwiseContr(design5))


trtVar = rowMeans(sapply(MS5[2,], function(x) sqrt(2*x/(6* contr.e))))


 trt.Est5 = apply(save.y5, 2, function(x) 2*tapply(x, design5$Trt, mean) %*% pairwiseContr(design5)* contr.e)

 p.value5 =  apply(trt.Est5, 2, function(x)  pt(-abs(x/trtVar), df = 5) * 2)

 power5 = apply(p.value5, 1, function(x) sum(x<0.05)/length(x))
 #power5 = apply(p.value5, 1, function(x) sum(x<0.05)/length(x))


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
#                 a  b  c  d  e  d
trt.Eff = rbind(c(3, 1, 3, 3, 5, 3),
                c(6, 1, 6, 6, 11, 6),
                c(3, 1, 3, 5, 3, 3),
                c(6, 1, 6, 11, 6, 6),
                c(1, 1, 1, 1, 1, 1),
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


 trt.Est6 =   2 *( trt.eff %*% pairwiseContr(design6))* contr.e

tPower6 = pt(qt(sig.level/2, lower.df, lower.tail = FALSE), lower.df, 
ncp =  abs(trt.Est6) /trtVar, lower.tail = FALSE)


      lower.df = 5
      lambda <- (groups - 1)  * 6  *  13160/15859 * (between.var/within.var)

fPower5 =  pf(qf(sig.level, groups - 1, lower.df , lower.tail = FALSE),
            groups - 1, lower.df , lambda, lower.tail = FALSE)


   contr.e = contrCompare(design5, pairwiseContr(design5))

trtVar = sqrt(2*(1 + 2*gamma.ani)/(6* contr.e))


 trt.Est5 =  2 *(trt.eff %*% pairwiseContr(design5))* contr.e

tPower5 =pt(qt(sig.level/2, lower.df, lower.tail = FALSE), lower.df, ncp = abs(trt.Est5) /trtVar, lower.tail = FALSE)

ref.temp6 = rbind(ref.temp6, round(c(i, gamma.ani, fPower6, tPower6), 3))

ref.temp5 = rbind(ref.temp5,  round(c(i, gamma.ani, fPower5, tPower5), 3))
}

}            
colnames(ref.temp6) = c("Trt.Eff", "gamma.ani", "F-test power",
  apply(states, 1, function(x) paste(x, sep = "", collapse = "-")))

colnames(ref.temp5) = c("Trt.Eff", "gamma.ani", "F-test power",
  apply(states, 1, function(x) paste(x, sep = "", collapse = "-")))

 eff = round(c(0, 0, 1/mean(1/contrCompare(design6, pairwiseContr(design6))), contrCompare(design6, pairwiseContr(design6))), 3)

ref.temp6 = rbind(eff, ref.temp6)


  eff =  round(c(0, 0, 1/mean(1/contrCompare(design5, pairwiseContr(design5))), contrCompare(design5, pairwiseContr(design5))), 3)

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


  contr6.e = contrCompare(design6, pairwiseContr(design6))
   contr5.e = contrCompare(design5, pairwiseContr(design5))

temp6 = data.frame()
temp5 = data.frame()

nSim = 10000

VC.resid = 1

gamma.run =  10
gamma.cag =  5

gamma.A = c(10^((-4:4)/2))
#                 a  b  c  d  e  d
trt.Eff = rbind(c(3, 1, 3, 3, 5, 3),
                c(6, 1, 6, 6, 11, 6),
                c(3, 1, 3, 5, 3, 3),
                c(6, 1, 6, 11, 6, 6),
                c(1, 1, 1, 1, 1, 1),
                c(1, 2, 3, 4, 5, 6),
                c(1, 5, 9, 13, 17, 21))

for(i in 1:nrow(trt.Eff)){
  cat("trt.Eff ", i, "\n")

  trt.eff = trt.Eff[i,]


    # create progress bar


test = function(j) {
    sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R",
    trace = FALSE)

    counter = 0

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



    trtVar = rowMeans(sapply(MS6[2, ], function(x) sqrt(2 * x/(6 * contr6.e))))


    trt.Est6 = apply(save.y6, 2, function(x) 2 * tapply(x, design6$Trt, mean) %*% pairwiseContr(design6) *
        contr6.e)

    p.value6 = apply(trt.Est6, 2, function(x) pt(-abs(x/trtVar), df = 6) * 2)

    power6 = apply(p.value6, 1, function(x) sum(p.adjust(x, "fdr") < 0.05)/length(x))



    trtVar = rowMeans(sapply(MS5[2, ], function(x) sqrt(2 * x/(6 * contr5.e))))


    trt.Est5 = apply(save.y5, 2, function(x) 2 * tapply(x, design5$Trt, mean) %*% pairwiseContr(design5) *
        contr5.e)

    p.value5 = apply(trt.Est5, 2, function(x) pt(-abs(x/trtVar), df = 5) * 2)

    power5 = apply(p.value5, 1, function(x) sum(p.adjust(x, "fdr") < 0.05)/length(x))
    # power5 = apply(p.value5, 1, function(x) sum(x<0.05)/length(x))


    return(c(gamma.A[j], sum(p.adjust(1 - pf(MS6[1, ]/MS6[2, ], 5, 6), "fdr") < 0.05)/nSim, power6,
        sum(p.adjust(1 - pf(MS5[1, ]/MS5[2, ], 5, 5), "fdr") < 0.05)/nSim, power5))
}

    sfExport("design.df", "gamma.A", "trt.eff", "VC.resid",  "design5",
    "design6", "sourceDir", "anovaTable6", "anovaTable5", "y6", "nSim",
    "pairwiseContr", "contr5.e", "contr6.e")

    temp = t(sfSapply(1:length(gamma.A),  fun = test))



    temp6 = rbind(temp6, cbind(i,temp[,1:17]))

    temp5 = rbind(temp5, cbind(i,temp[,c(1,18:33)]))

}
    
    
 sfStop()

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


