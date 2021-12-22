

orberg1 = data.frame( Run = factor(rep(1:6, each = 4)),
                      Tag = factor(rep(1:4, time = 6)),
                      Trt = factor( c(1,2,3,4,
                                      2,4,1,3,
                                      4,3,1,2,
                                      3,1,4,2,
                                      2,1,4,3,
                                      4,2,3,1)),
                      Sam = factor(rep(1:6, each = 4)))
                                      

summaryAovTwoPhase(orberg1,  blk.str2 = "Run", blk.str1 = "Trt:Sam",
trt.str = "Tag + Trt")


old = proc.time()
design.df = optCRD(nTrt = 4, bRep  = 6, tRep  = 1, nPlot = 4, iter  = 10000)
proc.time() - old

design1.df = optCRD(nTrt = 4, bRep  = 6, tRep  = 1, nPlot = 4, iter  = 1000)


design.summary.CRD(design1.df, simple = FALSE)


design.df = optCRD(nTrt = 4, bRep  = 8, tRep  = 1, nPlot = 8, iter  = 1000)

design.summary.CRD(design.df, simple = FALSE)


design1.df = optCRD(nTrt = 6, bRep  = 3, tRep  = 2, nPlot = 4, iter  = 1000)

design = design1.df

run.eff = rnorm(nlevels(design$Run), mean = 0, sd = sqrt(10))
ani.eff = rnorm(nlevels(design$Ani), mean = 0, sd = sqrt(5))
trt.eff = runif(nlevels(design$Trt), 0, 2)
tag.eff = runif(nlevels(design$Tag), 0, 1)
res.eff = rnorm(nrow(design), mean = 0, sd = sqrt(1))
gm = 10


y = gm + with(design, run.eff[Run] + ani.eff[Ani] +  tag.eff[Tag] + trt.eff[Trt]) + res.eff

aov.table = summaryAovTwoPhase(design.df = design, blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = y)

library(nlme)
m1.lme4 = lme(y ~ Tag + Trt,  random = ~ 1|Run + 1|Ani,
                       data = design)
 
summary(m1.lme4)

anova(m1.lme4)

















  
orberg1 = data.frame( Run = factor(rep(1:5, each = 4)),
                      Tag = factor(rep(1:4, time = 5)),
                      Trt = factor( c(1,2,3,4,
                                      2,3,4,1,
                                      3,4,1,2,
                                      4,1,2,3,
                                      2,1,4,3)))
                      
summaryAovOnePhase(design.df = orberg1, blk.str = "Run", 
                    trt.str = "Tag + Trt")
                                     
                                      
orberg1 = data.frame( Run = factor(rep(1:6, each = 4)),
                      Tag = factor(rep(1:4, time = 6)),
                      Trt = factor( c(1,2,3,4,
                                      2,3,4,1,
                                      3,4,1,2,
                                      4,1,2,3,
                                      2,1,4,3,
                                      4,2,3,1)))
                      
summaryAovOnePhase(design.df = orberg1, blk.str = "Run", 
                    trt.str = "Tag + Trt")
                      


orberg1 = data.frame( Run = factor(rep(1:3, each = 8)),
                      Tag = factor(rep(1:8, time = 3)),
                      Trt = factor( c(1,2,3,4,
                                      2,4,1,3,
                                      4,3,1,2,
                                      3,1,4,2,
                                      2,1,4,3,
                                      4,2,3,1)))
                                      
summaryAovOnePhase(design.df = orberg1, blk.str = "Run", 
                    trt.str = "Tag + Trt")







                                      