library(snowfall)   #Use parallel programming to speed up the computation by 4 times
 sfInit( parallel=TRUE, cpus=4 )

 (anovaTable14 = summaryAovTwoPhase(design.df = design14 ,blk.str2 = "Run",
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)

 (anovaTable16 = summaryAovTwoPhase(design.df = design16 ,blk.str2 = "Run",
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)


  contr16.e = contrCompare(design16, pairwiseContr(design16))
   contr14.e = contrCompare(design14, pairwiseContr(design14))

   trtEff16 = trtProj(design16)
   trtEff14 = trtProj(design14)

temp16 = data.frame()
temp14 = data.frame()

nSim = 10000

VC.resid = 1

gamma.run =  10
gamma.cag =  5

gamma.A = c(10^((-4:4)/2))
trt.Eff = rbind(c(3, 3, 3, 3, 1, 5, 3, 3),
                c(6, 6, 6, 6, 1, 11, 6, 6),
                c(3, 3, 1, 3, 7, 1, 3, 3),
                c(6, 6, 1, 6, 11, 6, 6, 6),
                c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                c(1, 2, 3,4,5,6, 7, 8))

for(i in 1:nrow(trt.Eff)){
  cat("trt.Eff ", i, "\n")

 
  trt.eff = trt.Eff[i,]

    # create progress bar


test = function(j) {
    #sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R",trace = FALSE)
     sourceDir(path = "C:/Users/Kevin/Dropbox/Packaging tools/packaging/infoDecompuTE/R",trace = FALSE)
    counter = 0

    # create progress bar
  
    save.y16 = rep(0, length(y16))
    save.y14 = rep(0, length(y16))

     MS14 = rep(0, 2)
            MS16 = rep(0, 2)
            #MSCRD = rep(0, nrow(anovaTableCRD))

             counter = 0
            # simulation


            while (counter < nSim) {


                run.eff = rnorm(nlevels(design.df$Run), mean = 0, sd = sqrt(10 * VC.resid))
                cag.eff = rnorm(nlevels(design.df$Cag), mean = 0, sd = sqrt(5 * VC.resid))
                ani.eff = rnorm(nlevels(design.df$Ani), mean = 0, sd = sqrt(gamma.A[j] * VC.resid))
                 tag.eff = runif(nlevels(design.df$Tag), 0, 1)
                res.eff = rnorm(nrow(design.df), mean = 0, sd = sqrt(VC.resid))
                gm = 10

                y14 = gm + with(design14, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] +
                        tag.eff[Tag] + trt.eff[Trt]) + res.eff

                y16 = gm + with(design16, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] +
                        tag.eff[Tag] + trt.eff[Trt]) + res.eff

              save.y16 = cbind(save.y16, y16)


              save.y14 = cbind(save.y14, y14)

                aov.table14 = summaryAovTwoPhase(design.df = design14, blk.str2 = "Run",
                blk.str1 = "Cag/Ani", trt.str = "Tag + Trt",
                  response = y14)

                aov.table16 = summaryAovTwoPhase(design.df = design16, blk.str2 = "Run",
                blk.str1 = "Cag/Ani", trt.str = "Tag + Trt",
                  response = y16)


                MS14 = cbind(MS14, as.numeric(aov.table14$A[10:11, "MS"]))
                MS16 = cbind(MS16, as.numeric(aov.table16$A[9:10, "MS"]))

                counter = counter + 1
            }


            MS14 = MS14[, -1]

            MS16 = MS16[, -1]

    save.y16 = save.y16[, -1]
    save.y14 = save.y14[, -1]



    trtVar = rowMeans(sapply(MS16[2, ], function(x) sqrt(2 * x/(8 * contr16.e))))


    trt.Est16 = apply(save.y16, 2, function(x) 2 * (trtEff16 %*% x) %*% pairwiseContr(design16) )

    p.value16 = apply(trt.Est16, 2, function(x) pt(-abs(x/trtVar), df = 16) * 2)

    power16 = apply(p.value16, 1, function(x) sum(p.adjust(x, "fdr") < 0.05)/length(x))



    trtVar = rowMeans(sapply(MS14[2, ], function(x) sqrt(2 * x/(8 * contr14.e))))


    trt.Est14 = apply(save.y14, 2, function(x) 2 *  (trtEff14 %*% x) %*% pairwiseContr(design14))

    p.value14 = apply(trt.Est14, 2, function(x) pt(-abs(x/trtVar), df = 14) * 2)

    power14 = apply(p.value14, 1, function(x) sum(p.adjust(x, "fdr") < 0.05)/length(x))


    return(c(gamma.A[j], sum(p.adjust(1 - pf(MS16[1, ]/MS16[2, ], 7, 16), "fdr") < 0.05)/nSim, power16,
        sum(p.adjust(1 - pf(MS14[1, ]/MS14[2, ], 7, 14), "fdr") < 0.05)/nSim, power14))
}

    sfExport("design.df", "gamma.A", "trt.eff", "VC.resid",  "design14",
    "design16", "sourceDir", "anovaTable16", "anovaTable14", "y16", "nSim",
    "pairwiseContr", "contr14.e", "contr16.e", "trtEff16", "trtEff14")

    temp = t(sfSapply(1:length(gamma.A),  fun = test))



    temp16 = rbind(temp16, cbind(i,temp[,1:30]))

    temp14 = rbind(temp14, cbind(i,temp[,c(1,31:59)]))

}


 sfStop()

 states = states[,1:2]
 
colnames(temp16) = c("Trt.Eff", "gamma.ani", "F-test power",
  apply(states, 1, function(x) paste(letters[x], sep = "", collapse = " - ")))

colnames(temp14) = c("Trt.Eff", "gamma.ani", "F-test power",
  apply(states, 1, function(x) paste(letters[x], sep = "", collapse = " - ")))

 eff = round(c(0, 0, 1/mean(1/contrCompare1(design16,  basicContr(design16)[,-8])),
  contrCompare(design16, pairwiseContr(design16))), 4)

temp16 = rbind(eff, temp16)


  eff =  round(c(0, 0, 1/mean(1/contrCompare1(design14,  basicContr(design14)[,-8])), 
  contrCompare(design14, pairwiseContr(design14))), 4)

temp14 = rbind(eff, temp14)


temp16[,"Trt.Eff"] = c(0,apply(trt.Eff[temp16[,"Trt.Eff"][-1],], 1, function(x) paste(x, sep = "", collapse = ",")))

temp14[,"Trt.Eff"] = c(0,apply(trt.Eff[temp14[,"Trt.Eff"][-1],], 1, function(x) paste(x, sep = "", collapse = ",")))

