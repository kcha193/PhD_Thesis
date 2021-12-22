
library(optimTE)



optRBD.combine = function(nTrt, bRep, nCag, tRep, nPlot, iter, confoundCag = FALSE) {
  
  design.df = optCRD(nTrt = nTrt, bRep = bRep, tRep = tRep, nPlot = nPlot,
                     iter = iter/10)
  n = nrow(design.df)
  nBlk = nlevels(design.df$Run)
  nPlot = nlevels(design.df$Tag)
  Run.mat = with(design.df, as.matrix(table(1:n, Run)))
  Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
  Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
  Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))
  Pb = projMat(Run.mat)
  Pb1 = projMat(Tag.mat)
  
  (ave.eff = test.CRD(X.trt = Trt.mat, (diag(n) - Pb) %*% (diag(n) - Pb1),
                      Rep = n/ncol(Trt.mat))$ave.eff)
  
  return(optRBD(nTrt = nTrt, bRep = bRep, nCag = nCag, tRep = tRep, nPlot = nPlot,
                iter = iter, upperValue = ave.eff, confoundCag = confoundCag))
} 



design.df = optRBD(nTrt = 3, bRep = 4, nCag = 2, tRep = 2, nPlot = 4, 
                           iter = 1000,
                           confoundCag = TRUE)

design1.df = optRBD(nTrt = 3, bRep = 4, nCag = 2, tRep = 2, nPlot = 4, 
                   iter = 100,
                   confoundCag = FALSE)

summaryAovTwoPhase(design1.df,
                   blk.str1 = "Cag/Ani", 
                   blk.str2 = "Run", 
                   trt.str = "Tag + Trt")


###############################################################################################

summaryAovTwoPhase(design.df,
                   blk.str1 = "Cag/Ani", 
                   blk.str2 = "Run", 
                   trt.str = "Tag + Trt")

bestDes <- design.df

################################################################################################

table(design.df$Cag, design.df$Ani)

table(design.df$Cag, design.df$Ani, design.df$Trt)


initialDes <- 
  cbind(design.df[,1:2], 
        design.df[order(design.df$Cag,
                        design.df$Ani,
                        design.df$Trt),-c(1:2)])

anovaTab <- 
summaryAovTwoPhase(initialDes,
                   blk.str1 = "Cag/Ani", 
                   blk.str2 = "Run", 
                   trt.str = "Tag + Trt")


resDF <- anovaTab$ANOVA[,1]
resDF <- resDF[trimws(names(resDF)) == "Residual"]

aveTrtEff <- anovaTab$Fixed$EF[, "eff.Trt"]
aveTrtEff <- aveTrtEff[trimws(names(aveTrtEff)) == "Trt"]

################################################################################################
#Generate 1,000,000 designs

library(parallel)

cl <- makeCluster(4)

n <- 1000000

old <- proc.time()

clusterExport(cl, "design.df")

design.list <- 
  parLapply(cl, 1:n, function(x){
         set.seed(x)
    library(infoDecompuTE)
          design <- cbind(design.df[,1:2], 
               design.df[sample(1:24),-c(1:2)])
         
          anovaTab <-
            summaryAovTwoPhase(
              design,
              blk.str1 = "Cag/Ani",
              blk.str2 = "Run",
              trt.str = "Tag + Trt", 
              decimal = TRUE, digits = 10
            )
          resDF <- anovaTab$ANOVA[, 1]
          resDF <- resDF[trimws(names(resDF)) == "Residual"]
          
          trtDF <- anovaTab$ANOVA[, 1]
          trtDF <- trtDF[trimws(names(trtDF)) == "Trt"]         
         
          aveTrtEff <- anovaTab$Fixed$EF[, "eff.Trt"]
          aveTrtEff <- aveTrtEff[trimws(names(aveTrtEff)) == "Trt"]
         
          list(design = design, 
               resDF = resDF, 
               trtDF=trtDF, 
               aveTrtEff = aveTrtEff)
         })

proc.time() - old

stopCluster(cl)


sample(1:24)


#save(design.list, bestDes, file = "designs.Rdata")

summaryAovTwoPhase(
  bestDes,
  blk.str1 = "Cag/Ani",
  blk.str2 = "Run",
  trt.str = "Tag + Trt"
)



summaryAovTwoPhase(
  design.list[[8]]$design,
  blk.str1 = "Cag/Ani",
  blk.str2 = "Run",
  trt.str = "Tag + Trt"
)

load("designs.Rdata")

trtDFs <- as.numeric(sapply(design.list, function(x) rev(x$trtDF)[1]) )
range(trtDFs)

resDFs <- as.numeric(sapply(design.list, function(x) rev(x$resDF)[2]) )
range(resDFs)

which(resDFs == max(resDFs))

design.summary.RBD(design.list[[56648]]$design)

design.summary.RBD(design.list[[141791]]$design)

aveTrtEffs <- as.numeric(sapply(design.list, function(x) rev(x$aveTrtEff)[1]) )
range(aveTrtEffs)
max(aveTrtEffs)

which(aveTrtEffs == max(aveTrtEffs))

design.summary.RBD(design.list[[25480]]$design)
design.summary.RBD(design.list[[38650]]$design)
design.summary.RBD(design.list[[47237]]$design)
design.summary.RBD(design.list[[54865]]$design)
design.summary.RBD(design.list[[77264]]$design)

design.summary.RBD(bestDes)

bestresDFs <-  max(resDFs)
besAveTrtEffs <-  max(aveTrtEffs)

rm(design.list, aveTrtEffs, trtDFs, resDFs)
gc()


##################################################

library(infoDecompuTE)

design.df <- bestDes


table(design.df$Cag, design.df$Ani)

table(design.df$Cag, design.df$Ani, design.df$Trt)


design.listNew <- list()

pb <- txtProgressBar(min = 2e6, max = 3e6, style = 3)

for(i in 2e6:3e6){

  # update progress bar
  setTxtProgressBar(pb, i)
  
  set.seed(i)
  design <- cbind(design.df[,1:2], 
                  design.df[sample(1:24),-c(1:2)])
  
  anovaTab <-
    summaryAovTwoPhase(
      design,
      blk.str1 = "Cag/Ani",
      blk.str2 = "Run",
      trt.str = "Tag + Trt", 
      decimal = TRUE, digits = 10
    )
  
  
  resDF <- anovaTab$ANOVA[, 1]
  resDF <- resDF[trimws(names(resDF)) == "Residual"]

  aveTrtEff <- anovaTab$Fixed$EF[, "eff.Trt"]
  aveTrtEff <- aveTrtEff[trimws(names(aveTrtEff)) == "Trt"]
  
  
  if(bestresDFs == rev(resDF)[2] || besAveTrtEffs == rev(aveTrtEff)[1]) {
    
    cat("\n Saving this design", i, "\n")
  
    design.listNew <- list(design.listNew,
                        
      list(design = design, 
           resDF = resDF, 
           trtDF=trtDF, 
           aveTrtEff = aveTrtEff)
    )
  }
  
}
close(pb)

design.listNew

#########################################################################################################

Tray_Tag <- 
expand.grid(Tag1= c(1,2),
            Tag2= c(1,2),
            Tag3= c(1,2),
            Tag4= c(1,2))

Tray_Tag <- Tray_Tag[-c(1,nrow(Tray_Tag)),]


Tray_Tag <- 
  Tray_Tag[apply(Tray_Tag, 1, 
               function(x) 
                 table(x)[1] == table(x)[2]  ),]


(Tray_Tag <- Tray_Tag[c(3,5,6),])


Tray_Run <- 
expand.grid(Run1= c(1,2),
            Run2= c(1,2),
            Run3= c(1,2),
            Run4= c(1,2),
            Run5= c(1,2),
            Run6= c(1,2))

Tray_Run <- Tray_Run[-c(1,nrow(Tray_Run)),]

Tray_Run <- 
  Tray_Run[apply(Tray_Run, 1, 
                 function(x) 
                   table(x)[1] == table(x)[2]  ),]


(Tray_Run <- Tray_Run[Tray_Run[,1]==1,])



###############################################################################

Trt_Tray <- 
  expand.grid(Sam1=  c(1,2,3),
              Sam2=  c(1,2,3),
              Sam3=  c(1,2,3),
              Sam4=  c(1,2,3),
              Sam5=  c(1,2,3),
              Sam6=  c(1,2,3),
              Sam7=  c(1,2,3),
              Sam8=  c(1,2,3),
              Sam9=  c(1,2,3),
              Sam10= c(1,2,3),
              Sam11= c(1,2,3),
              Sam12= c(1,2,3))

index <- apply(Trt_Tray, 1, function(x) length(table(x)) == 3)

Trt_Tray <-   Trt_Tray[index,]

index2 <- apply(Trt_Tray, 1, 
                function(x) 
                  table(x)[1] == table(x)[2] & 
                  table(x)[1] ==table(x)[3] & 
                  table(x)[2] == table(x)[3])

Trt_Tray <-   Trt_Tray[index2,]


Trt_Tray <- 
  Trt_Tray[order(Trt_Tray[,1], Trt_Tray[,2], Trt_Tray[,3], 
               Trt_Tray[,4], Trt_Tray[,5], Trt_Tray[,6], 
               Trt_Tray[,7], Trt_Tray[,8], Trt_Tray[,9], 
               Trt_Tray[,10], Trt_Tray[,11], Trt_Tray[,12]), ]

dim(Trt_Tray)


###############################################################################

Plant_Trt <- 
  expand.grid(Sam1= c(1,4),
              Sam2= c(1,4),
              Sam3= c(1,4),
              Sam4= c(1,4))

Plant_Trt <- Plant_Trt[-c(1,nrow(Plant_Trt)),]


Plant_Trt <- 
  Plant_Trt[apply(Plant_Trt, 1, 
                 function(x) 
                   table(x)[1] == table(x)[2]  ),]


(Plant_Trt <- Plant_Trt[c(3,5,6),])


############################################

temp <- rep(1:34650, 34650)
temp1 <- rep(1:34650, each = 34650)




########################################################################
displayDes <- function(des){
  
  des <- matrix(paste0(des$Tray, des$Plant, des$Trt), ncol = 4, nrow = 6, byrow = TRUE)
  rownames(des) <- 1:6
  colnames(des) <- 114:117
  
  knitr::kable(des,  row.names = TRUE)
}

##############################################################################
library(infoDecompuTE)
initalDes <- 
data.frame(Run = rep(1:6, each = 4), 
           Tag = rep(114:117, 6), 
           Tray = rep(1:2, each = 12), 
           Plant = rep(LETTERS[1:12], each= 2), 
           Trt = rep(letters[1:3], each = 2))

summaryAovTwoPhase(
  initalDes,
  blk.str1 = "Tray/Plant",
  blk.str2 = "Run",
  trt.str = "Tag + Trt"
)

############################################################################


# Tray by Tags

initalDes1 <- 
  cbind(initalDes[,1:2], initalDes[as.numeric(t(matrix(1:24, ncol = 2))) , 3:5])


c(1,2,3,4) # 1212
c(1,3,2,4) # 1122
c(1,3,4,2) # 1221



initalDes2 <- 
  cbind(initalDes1[,1:2], 
        initalDes1[as.numeric(t(matrix(1:24, ncol = 4, byrow = TRUE)[,c(1,2,3,4)])) , 3:5])

Trt_Tray_index <- 1

initalDes2[initalDes2$Tray == 1,"Trt"] <- 
  letters[as.numeric(Trt_Tray[Trt_Tray_index,])]
initalDes2[initalDes2$Tray == Trt_Tray_index,"Trt"] <- 
  letters[as.numeric(Trt_Tray[Trt_Tray_index,])]


Plant_Trt_index <- 1

initalDes2[initalDes2$Tray == 1 & initalDes2$Trt == "a","Plant"] <- 
  LETTERS[as.numeric(Plant_Trt[Plant_Trt_index,])]
initalDes2[initalDes2$Tray == 1 & initalDes2$Trt == "b","Plant"] <- 
  LETTERS[as.numeric(Plant_Trt[Plant_Trt_index,]) + 1]
initalDes2[initalDes2$Tray == 1 & initalDes2$Trt == "c","Plant"] <-  
  LETTERS[as.numeric(Plant_Trt[Plant_Trt_index,]) + 2]
initalDes2[initalDes2$Tray == 2 & initalDes2$Trt == "a","Plant"] <- 
  LETTERS[as.numeric(Plant_Trt[Plant_Trt_index,]) + 3]
initalDes2[initalDes2$Tray == 3 & initalDes2$Trt == "b","Plant"] <-  
  LETTERS[as.numeric(Plant_Trt[Plant_Trt_index,]) + 4]
initalDes2[initalDes2$Tray == 4 & initalDes2$Trt == "c","Plant"] <- 
  LETTERS[as.numeric(Plant_Trt[Plant_Trt_index,]) + 5]


displayDes(initalDes2)

summaryAovTwoPhase(
  initalDes2,
  blk.str1 = "Tray/Plant",
  blk.str2 = "Run",
  trt.str = "Tag + Trt"
)



displayDes(initalDes1)



summaryAovTwoPhase(
  initalDes1,
  blk.str1 = "Tray/Plant",
  blk.str2 = "Run",
  trt.str = "Tag + Trt"
)









# Tray by Runs

as.numeric(t(matrix(1:24, ncol = 4)))




initalDes[ , 3:5]

c(1,4,5,6,2,3) #122211
c(1,4,5,2,6,3) #122121
c(1,4,2,5,6,3) #121221
c(1,2,4,5,6,3) #112221
c(1,4,5,2,3,6) #122112
c(1,4,2,5,3,6) #121212
c(1,2,4,5,3,6) #112212
c(1,4,2,3,5,6) #121122
c(1,2,4,3,5,6) #112122
c(1,2,3,4,5,6) #111222


des1 <- 
  initalDes[as.numeric(t(matrix(1:24, ncol = 4, byrow = TRUE)[c(1,4,5,6,2,3) ,])) , 3:5]





#####################################


temp <- 
  expand.grid(Sam1=  LETTERS[c(1,2,3,4)],
              Sam2=  LETTERS[c(1,2,3,4)],
              Sam3=  LETTERS[c(1,2,3,4)],
              Sam4=  LETTERS[c(1,2,3,4)],
              Sam5= LETTERS[c(1,2,3,4)],
              Sam6=  LETTERS[c(1,2,3,4)],
              Sam7=  LETTERS[c(1,2,3,4)],
              Sam8=  LETTERS[c(1,2,3,4)])


index <- apply(temp, 1, function(x) length(table(x)) == 4)


temp <-   temp[index,]

index2 <- apply(temp, 1, 
                function(x) 
                  table(x)[1] == table(x)[2] & 
                  table(x)[2] ==table(x)[3] & 
                  table(x)[3] == table(x)[4])

temp <- temp[index2,]

dim(temp)



temp <- 
  expand.grid(Sam1=  LETTERS[c(1,2,3)],
              Sam2=  LETTERS[c(1,2,3)],
              Sam3=  LETTERS[c(1,2,3)],
              Sam4=  LETTERS[c(1,2,3)],
              Sam5= LETTERS[c(1,2,3)],
              Sam6=  LETTERS[c(1,2,3)])


index <- apply(temp, 1, function(x) length(table(x)) == 3)


temp <-   temp[index,]

index2 <- apply(temp, 1, 
                function(x) 
                  table(x)[1] == table(x)[2] & 
                  table(x)[1] ==table(x)[3] & 
                  table(x)[2] == table(x)[3])

temp <- temp[index2,]

dim(temp)





temp1 <- 
  expand.grid(Sam1=  LETTERS[c(1,2,3)],
              Sam2=  LETTERS[c(1,2,3)],
              Sam3=  LETTERS[c(1,2,3)],
              Sam4=  LETTERS[c(1,2,3)],
              Sam5= LETTERS[c(1,2,3)],
              Sam6=  LETTERS[c(1,2,3)],
              Sam7=  LETTERS[c(1,2,3)],
              Sam8=  LETTERS[c(1,2,3)],
              Sam8=  LETTERS[c(1,2,3)])


index <- apply(temp1, 1, function(x) length(table(x)) == 3)


temp1 <-   temp1[index,]

index2 <- apply(temp1, 1, 
                function(x) 
                  table(x)[1] == table(x)[2] & 
                  table(x)[1] ==table(x)[3] & 
                  table(x)[2] == table(x)[3])

temp1 <- temp1[index2,]

dim(temp1)




############################################################





























