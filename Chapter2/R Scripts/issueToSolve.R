

des1 <- local({
  
  Run <- factor(rep(1:6, each=4))
  Animal <- factor(c(1,6,5,4,3,2,
                     2,1,6,5,4,3,
                     3,2,1,6,5,4,
                     4,3,2,1,6,5))
  Tag <- factor(rep(114:117, times=6))
  Trt <- factor(rep(LETTERS[1:2], times=3)[Animal])
  Time <- factor(rep(c(1:3,3), each=6))
  Sample <- factor(rep(c(1:3,3), each=6))
  Subsample <- factor(rep(c(1,1,1,2), times=6))
  
  data.frame(Run, Tag, Animal, Trt, Sample, Time, Subsample)
  
}) 
des1


library(infoDecompuTE)

# ANOVA for Phase 1 design
des1Ph1 <- with(des1, des1[,c("Animal", "Sample", "Subsample",
                              "Tag", "Trt", "Time")])
des1Ph1$obs <- 1:nrow(des1Ph1)

with(des1Ph1, des1Ph1[order(Animal, Time, Subsample),])
summaryAovOnePhase(design.df=des1Ph1,
                   blk.str = "Animal/Sample", 
                   trt.str="Trt*Time")


# $Fixed$EF
# eff.Trt eff.Time eff.Trt$Time
# Between Animal                              
# Trt         1                1/9         
# Within Animal                               
# Time                1                    
# Trt$Time                     1           




des1$Inter <-interaction(des1$Trt, des1$Time, sep = "")

summaryAovOnePhase(design.df=des1Ph1,
                   blk.str = "Animal/Sample", 
                   trt.str="Trt + Time + Inter")


summaryAovTwoPhase(design.df=des1, blk.str2 = "Run",
                   blk.str1 = "Animal/Sample", 
                   trt.str="Tag +Trt *Time", decimal = TRUE)

des1$Inter <-interaction(des1$Trt, des1$Time, sep = "")

summaryAovTwoPhase(design.df=des1, blk.str2 = "Run",
                   blk.str1 = "Animal/Sample", 
                   trt.str="Tag + Trt + Time + Inter", decimal = TRUE)


library(dae)
##obtain projectors using projs.structure
Q.unit <- projs.structure(~ Animal/Sample, data = des1Ph1)
Q.trt <- projs.structure(~ Trt*Time, data = des1Ph1)

proj2.efficiency(Q.unit[["Animal"]], Q.trt[["Trt"]])
proj2.efficiency(Q.unit[["Animal"]], Q.trt[["Time"]])
proj2.efficiency(Q.unit[["Animal"]], Q.trt[["Trt:Time"]])

proj2.efficiency(Q.unit[["Animal:Sample"]], Q.trt[["Trt"]])
proj2.efficiency(Q.unit[["Animal:Sample"]], Q.trt[["Time"]])
proj2.efficiency(Q.unit[["Animal:Sample"]], Q.trt[["Trt:Time"]])

##obtain combined decomposition and summarize
unit.trt.p2canon <- projs.2canon(Q.unit, Q.trt)
summary(unit.trt.p2canon)


des1$Inter <-interaction(des1$Trt, des1$Time, sep = "")

summaryAovOnePhase(design.df=des1Ph1,
                   blk.str = "Animal/Sample", 
                   trt.str="Trt + Time + Inter")


##obtain projectors using projs.structure
Q.unit <- projs.structure(~ Animal/Sample, data = des1Ph1)
Q.trt <- projs.structure(~ Trt + Time + Inter, data = des1Ph1)

proj2.efficiency(Q.unit[["Animal"]], Q.trt[["Trt"]])
proj2.efficiency(Q.unit[["Animal"]], Q.trt[["Time"]])
proj2.efficiency(Q.unit[["Animal"]], Q.trt[["Inter"]])

proj2.efficiency(Q.unit[["Animal:Sample"]], Q.trt[["Trt"]])
proj2.efficiency(Q.unit[["Animal:Sample"]], Q.trt[["Time"]])
proj2.efficiency(Q.unit[["Animal:Sample"]], Q.trt[["Inter"]])

##obtain combined decomposition and summarize
unit.trt.p2canon <- projs.2canon(Q.unit, Q.trt)
summary(unit.trt.p2canon)




summaryAovOnePhase(design.df=des1Ph1, blk.str = "obs", 
                   trt.str="Trt")

summaryAovOnePhase(design.df=des1Ph1, blk.str = "Animal", 
                   trt.str="Trt*Time")


table(des1Ph1$Trt, des1Ph1$Tag)


inter <-interaction(des1Ph1$Trt, des1Ph1$Time)

interDes <- table(1:length(inter), inter)
trtDes <- table(1:length(inter), des1Ph1$Trt)
timeDes <- table(1:length(inter), des1Ph1$Time)

aniDes <- table(1:length(inter), des1Ph1$Animal)
samDes <- table(1:length(inter), des1Ph1$Sample)



inforMat <- t(trtDes) %*%  (projMat(aniDes) - K(24))  %*% trtDes

inforMat <- t(interDes) %*%  (projMat(aniDes) - K(24))  %*% interDes

inforMat <- t(interDes) %*%  (identityMat(24) - projMat(aniDes))  %*% interDes




inforMat <- T$Trt %*% t(interDes) %*%  (projMat(aniDes) - K(24))  %*% interDes %*%  T$Trt
r <- diag(1/sqrt(rep(12, 6)))


eigen(inforMat)

inforMat <- t(interDes) %*%  (projMat(aniDes) - K(24))  %*% interDes

eigen(inforMat)


r <- diag(1/sqrt(as.numeric(table(inter))))

eigen(r %*% inforMat %*% r)


inforMat <-T$`Trt:Time` %*%  t(interDes) %*% projMat(samDes)%*% (identityMat(24) -projMat(aniDes)) %*% 
  projMat(samDes) %*% interDes %*%  T$`Trt:Time`


r <- diag(1/sqrt(as.numeric(table(inter))))

eigen(r %*% inforMat %*% r)$va








