

# Load in the designs #####

load("designTag2.Rdata")

desTag <- designTag1.list[[854022]]$design

rm(designTag1.list)


load("designRun1.Rdata")

desRun <- designRun1.list[[1357108]]$design

rm(designRun1.list)

load("designRunTag.Rdata")

desRunTag <- designRunTag.list[[673808]]$design

rm(designRunTag.list)


library(infoDecompuTE)
# Sim studies



getVarTrt <- function(des){

  des$Run <- factor(des$Run)
  des$Tag <- factor(des$Tag)
  des$Tray <- factor(des$Tray)
  des$Plant <- factor(des$Plant)
  des$Trt <- factor(des$Trt)
  
  n <- nrow(des)
  sigmaS <- 1
  sigmaP <- 5
  sigmaT <- 10
  sigmaR <- 100
  
  xTag <- c(1, 1, 1, 1)
  xTrt <- c(1, 6, 10)
  
  set.seed(12)
  y <- sapply(1:1000, function(x) {
    1000 + rnorm(nlevels(des$Run), sd = sigmaR)[des$Run] +
      rnorm(nlevels(des$Tray), sd = sigmaT)[des$Tray] +
      rnorm(nlevels(des$Plant), sd = sigmaP)[des$Plant] +
      rnorm(n, sd = sigmaS) +
      xTag[des$Tag] + xTrt[des$Trt]
  })
  
  aov.fit <- apply(y, 2, function(x)
    aov(x ~ Tag + Trt + Error(Run + Tray / Plant), des))
  
  aov.fit
}


VarTrtTag <- getVarTrt(desTag)


summaryAovTwoPhase(
  desTag,
  blk.str1 = "Tray/Plant",
  blk.str2 = "Run",
  trt.str = "Tag + Trt"
)

VarTrtTagEst <- 
  sapply(VarTrtTag, function(x) 
    rev(summary(x$`Tray:Plant`)[[1]][,3])[1]/8/(15/16))

VarTrtTagEstDiff <- 
  sapply(VarTrtTag, function(x)
    summary.lm(x$`Tray:Plant`)$coef[,1])

VarTrtTagEstSE <- 
  sapply(VarTrtTag, function(x)
    summary.lm(x$`Tray:Plant`)$coef[,2])



summary(VarTrtTagEst)

VarTrtRun <- getVarTrt(desRun)

summaryAovTwoPhase(
  desRun,
  blk.str1 = "Tray/Plant",
  blk.str2 = "Run",
  trt.str = "Tag + Trt"
)


VarTrtRunEst <- 
  sapply(VarTrtRun, function(x) 
    rev(summary(x$`Tray:Plant`)[[1]][,3])[1]/8/(195/224))

VarTrtRunEstDiff <- 
  sapply(VarTrtRun, function(x)
    summary.lm(x$`Tray:Plant`)$coef[2:3,1])

VarTrtRunEstSE <- 
  sapply(VarTrtRun, function(x)
    summary.lm(x$`Tray:Plant`)$coef[2:3,2])

summary(VarTrtRunEst)


VarTrtRunTag <- getVarTrt(desRunTag)

summaryAovTwoPhase(
  desRunTag,
  blk.str1 = "Tray/Plant",
  blk.str2 = "Run",
  trt.str = "Tag + Trt"
)

VarTrtRunTagEst <- 
  sapply(VarTrtRunTag,function(x) 
    rev(summary(x$`Tray:Plant`)[[1]][,3])[1]/8/( 3/4))

summary(VarTrtRunTagEst)

VarTrtRunTagEstDiff <- 
  sapply(VarTrtRunTag, function(x)
    summary.lm(x$`Tray:Plant`)$coef[2:3,1])

VarTrtRunTagEstSE <- 
  sapply(VarTrtRun, function(x)
    summary.lm(x$`Tray:Plant`)$coef[2:3,2])


# Estimates Mean and SD #####

# Trt B - A = 5  
# Trt B - A = 9  

apply(VarTrtTagEstDiff, 1, mean)

VarTrtTagEstBias <- VarTrtTagEstDiff
VarTrtTagEstBias[1,] <- VarTrtTagEstBias[1,]- 5
VarTrtTagEstBias[2,] <- VarTrtTagEstBias[2,]- 9

apply(VarTrtTagEstBias, 1, mean)


apply(VarTrtTagEstDiff, 1, sd)

apply(VarTrtRunEstDiff, 1, mean)

VarTrtRunEstBias <- VarTrtRunEstDiff
VarTrtRunEstBias[2,] <- VarTrtRunEstBias[2,]- 5
VarTrtRunEstBias[3,] <- VarTrtRunEstBias[3,]- 9

apply(VarTrtRunEstBias, 1, mean)

apply(VarTrtRunEstDiff, 1, sd)

apply(VarTrtRunTagEstDiff, 1, mean)

VarTrtRunTagEstBias <- VarTrtRunTagEstDiff
VarTrtRunTagEstBias[2,] <- VarTrtRunTagEstBias[2,]- 5
VarTrtRunTagEstBias[3,] <- VarTrtRunTagEstBias[3,]- 9

apply(VarTrtRunTagEstBias, 1, mean)

apply(VarTrtRunTagEstDiff, 1, sd)



apply(VarTrtTagEstSE, 1, mean)
apply(VarTrtRunEstSE, 1, mean)
apply(VarTrtRunTagEstSE, 1, mean)



boxplot(VarTrtTagEstSE[1,])



#Disaply the result ####




library(tidyverse)

data.frame(Tag = VarTrtTagEst, 
           Run = VarTrtRunEst, 
           RunTag = VarTrtRunTagEst) %>% summary() %>% 
  knitr::kable()


data.frame(Tag = VarTrtTagEst, 
           Run = VarTrtRunEst, 
           RunTag = VarTrtRunTagEst) %>% 
  gather(key = DesignType, val = Var) %>% 
  mutate(DesignType= factor(DesignType, 
                            levels = c("Tag", "Run", "RunTag"))) %>% 
  ggplot(aes(x = DesignType, y = Var)) + 
  geom_boxplot() + theme_bw()







