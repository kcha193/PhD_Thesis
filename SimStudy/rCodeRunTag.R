library(infoDecompuTE)
library(MASS)
library(optimTE)

initalDes <- 
  data.frame(Run = rep(1:6, each = 4), 
             Tag = rep(114:117, 6), 
             Tray = rep(1:2, each = 12), 
             Plant = rep(LETTERS[1:12], each= 2), 
             Trt = rep(letters[1:3], each = 2))

initalDes1 <- 
  cbind(initalDes[,1:2], initalDes[as.numeric(t(matrix(1:24, ncol = 2))) , 3:5])


initalDes2 <- 
  cbind(initalDes1[,1:2], 
        initalDes1[as.numeric(t(matrix(1:24, ncol = 4, byrow = TRUE)[,c(1,2,3,4)]))[
          c(1:4, 8:5, 9:12, 16:13, 17:20, 24:21)] , 3:5])


displayDes(initalDes2)


##############################################################################
design.df <- initalDes2

n = nrow(design.df)
nBlk = 6
nPlot = 4
nCag = 2

nAni = 6
Run.mat = with(design.df, as.matrix(table(1:n, Run)))
Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))

C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))


C.plot = (identityMat(nPlot) - K(nPlot)) %x% K(3)
C.trt = K(nPlot) %x% (identityMat(3) - K(3))



Pb = projMat(Run.mat)
Pb1 = projMat(Tag.mat)

blk.proj = (identityMat(n) - Pb) #%*% (identityMat(n) - Pb1)




##########################################################################################

library(parallel)

old <- proc.time()

cl <- makeCluster(4)

clusterExport(cl, c( "index", "initalDes2", "Trt_Tray" ,"Plant_Trt", 
                     "blk.proj", "matMulti1", "C.ani", "C.cage", 
                     "matMulti", "Tag.mat", "Pb1", "C.trt", "C.plot"))

iteration <- 2*10^6

designRunTag.list <-
  parLapply(cl, 1:iteration, function(x) {
    
    library(infoDecompuTE)
    library(MASS)
    library(optimTE)
    
    index <- sample(1:34650, 2)
    
    Plant_Trt_index <- sample(1:3, size = 1)
    
    initalDes2[initalDes2$Tray == 1, "Trt"] <-
      letters[as.numeric(Trt_Tray[index[1], ])]
    # letters[as.numeric(Trt_Tray[4488, ])]
    
    initalDes2[initalDes2$Tray == 2, "Trt"] <-
      letters[as.numeric(Trt_Tray[index[2], ])]
    
    initalDes2[initalDes2$Tray == 1 & initalDes2$Trt == "a", "Plant"] <-
      LETTERS[as.numeric(Plant_Trt[Plant_Trt_index, ])]
    
    initalDes2[initalDes2$Tray == 1 & initalDes2$Trt == "b", "Plant"] <-
      LETTERS[as.numeric(Plant_Trt[Plant_Trt_index, ]) + 1]
    
    initalDes2[initalDes2$Tray == 1 &
                 initalDes2$Trt == "c", "Plant"] <-
      LETTERS[as.numeric(Plant_Trt[Plant_Trt_index, ]) + 2]
    
    initalDes2[initalDes2$Tray == 2 & initalDes2$Trt == "a", "Plant"] <-
      LETTERS[as.numeric(Plant_Trt[Plant_Trt_index, ]) + 6]
    
    initalDes2[initalDes2$Tray == 2 &
                 initalDes2$Trt == "b", "Plant"] <-
      LETTERS[as.numeric(Plant_Trt[Plant_Trt_index, ]) + 7]
    
    initalDes2[initalDes2$Tray == 2 & initalDes2$Trt == "c", "Plant"] <-
      LETTERS[as.numeric(Plant_Trt[Plant_Trt_index, ]) + 8]
    
    design <- initalDes2
    
    Ani.mat = with(design, as.matrix(table(1:24, paste(Tray, Plant, sep = ""))))
    #Trt.mat = with(design, as.matrix(table(1:24, Trt)))
    Trt.mat = with(design, as.matrix(table(1:24, interaction(Trt, Tag,  sep = ""))))
    
    blk.projCage <- blk.proj - 
      matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)), C.cage, Ani.mat)
    
    avePlantEff <-  test.RBD(
      X.trt =Ani.mat,
      blk.proj = blk.projCage,
      C.trt.mat = C.ani,
      Rep = 2
    )$ave.eff
    
    PP1 <- matMulti1(blk.projCage, ginv(matMulti(blk.projCage, Ani.mat, C.ani)), C.ani, Ani.mat)
    
    
    PP_tag <- matMulti1(PP1, ginv(matMulti(PP1, Trt.mat, C.plot)), C.plot, Trt.mat)
    
    PP_trt <- matMulti1(PP1 - PP_tag, 
                        ginv(matMulti(PP1 - PP_tag, Trt.mat, C.trt)), C.trt, Trt.mat)
    
    
    aveEffDf <- as.numeric(test.RBD(X.trt = Trt.mat, 
                                    PP_trt, C.trt.mat = C.trt, Rep = 2)[c("ave.eff", "nCan")])
    
    
    resDF <- tr(PP1) - aveEffDf[2] -  test.RBD(X.trt = Tag.mat, PP1, Rep = 8)$nCan
    
    
    list(
      design = design,
      avePlantEff = avePlantEff,
      resDF = resDF,
      trtDF = aveEffDf[2],
      aveTrtEff = aveEffDf[1]
    )
  })

proc.time() - old

stopCluster(cl)


save(designRunTag.list,  file = "designRunTag.Rdata")

load("designRunTag.Rdata")

avePlantEffs <- as.numeric(sapply(designRunTag.list, function(x) x$avePlantEff ))

trtDFs <- as.numeric(sapply(designRunTag.list, function(x) x$trtDF) )

resDFs <- as.numeric(sapply(designRunTag.list, function(x) x$resDF))

aveTrtEffs <- as.numeric(sapply(designRunTag.list, function(x) x$aveTrtEff) )

sum(avePlantEffs == 1)

sum(avePlantEffs == 1)/2000000

table(trtDFs[avePlantEffs == 1])

table(round(resDFs))

table(round(resDFs[avePlantEffs == 1 & trtDFs == 2]))


which(avePlantEffs == 1 & trtDFs == 2 & resDFs >= 5)

aveTrtEffs[which(avePlantEffs == 1 & trtDFs == 2 & resDFs >= 5)]


summaryAovTwoPhase(
  designRunTag.list[[673808]]$design,
  blk.str1 = "Tray/Plant",
  blk.str2 = "Run",
  trt.str = "Tag + Trt", latex = TRUE, decimal = TRUE
)

table(designRunTag.list[[673808]]$design$Tag, designRunTag.list[[673808]]$design$Trt)

displayDes(designRunTag.list[[673808]]$design)

rm(designRunTag.list)
gc()



summaryAovTwoPhase(
  designRunTag.list[[1505960]]$design,
  blk.str1 = "Tray/Plant",
  blk.str2 = "Run",
  trt.str = "Tag + Trt"
)













