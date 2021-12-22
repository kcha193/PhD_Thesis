


library(xtable)
library(optimTE)



des <- read.csv("OptimalDesigns/CRD/designCRD8824.csv")

summaryAovTwoPhase(
  design.df = des,
  blk.str2 = "Run",
  blk.str1 = "Ani",
  trt.str = "Tag + Trt",
  decimal = TRUE,
  digits = 4
)

des <- read.csv("OptimalDesigns/CRD/designCRD8228.csv")

summaryAovTwoPhase(
  design.df = des,
  blk.str2 = "Run",
  blk.str1 = "Ani",
  trt.str = "Tag + Trt",
  decimal = TRUE,
  digits = 4
)


out <- matrix(paste0(des$Cag), ncol = nlevels(factor(des$Tag)), 
              nrow = max(des$Run), byrow = TRUE)
out

summaryAovOnePhase(
  design.df = des,
  blk.str = "Cag/Ani",
  trt.str = "Trt",
  decimal = TRUE,
  digits = 4
)$A

summaryAovTwoPhase(
  design.df = des,
  blk.str2 = "Run",
  blk.str1 = "Cag/Ani",
  trt.str = "Tag + Trt",
  decimal = TRUE,
  digits = 4
)






des <- read.csv("OptimalDesigns/RCBD/designRCBD22224.csv")


des$ANI <- LETTERS[interaction(des$Ani, des$Cag)]


out <- matrix(paste0(des$Cag, des$ANI, des$Trt), ncol = 4, 
              nrow = 2, byrow = TRUE)

xtable(out)


des <- read.csv("OptimalDesigns/RCBD/designRCBD63324.csv")

interaction(des$Ani, des$Cag)

(des$ANI <- LETTERS[interaction(des$Ani, des$Cag)])

sort(unique(paste0(des$Cag, des$ANI)))

out <-
  matrix(
    paste0(des$Cag, des$ANI, des$Trt),
    ncol = nlevels(factor(des$Tag)),
    nrow = max(des$Run),
    byrow = TRUE
  )

out

xtable(out)

summaryAovOnePhase(
  design.df = des,
  blk.str = "Cag/ANI",
  trt.str = "Trt",
  decimal = TRUE,
  digits = 4
)

summaryAovTwoPhase(
  design.df = des,
  blk.str2 = "Run",
  blk.str1 = "Cag/ANI",
  trt.str = "Tag + Trt",
  decimal = TRUE,
  digits = 4
)



multiLetters(1:30)

list.files("OptimalDesigns/RCBD/")

read.csv("OptimalDesigns/RCBD/")


des <- read.csv("OptimalDesigns/BIBD/designBIBD65624.csv")

interaction(des$Ani, des$Cag)

(des$ANI <-multiLetters(as.numeric(interaction(des$Ani, des$Cag))))

sort(unique(paste0(des$Cag, des$ANI)))

out1 <- matrix(sort(unique(paste0(des$ANI, des$Trt))), 
              nrow = max(des$Cag), byrow = TRUE)


xtable(out1)

out <- matrix(paste0(des$Cag, des$ANI, des$Trt), ncol = nlevels(factor(des$Tag)), 
              nrow = max(des$Run), byrow = TRUE)

xtable(out)

summaryAovOnePhase(
  design.df = des,
  blk.str = "Cag/ANI",
  trt.str = "Trt",
  decimal = TRUE,
  digits = 4
)$A


summaryAovOnePhase(
  design.df = des,
  blk.str = "ANI",
  trt.str = "Trt",
  decimal = TRUE,
  digits = 4
)$A

summaryAovOnePhase(
  design.df = des,
  blk.str = "Run",
  trt.str = "Cag/ANI",
  decimal = TRUE,
  digits = 4
)

summaryAovOnePhase(
  design.df = des,
  blk.str = "Run",
  trt.str = "Cag/Ani"
)

summaryAovOnePhase(
  design.df = des,
  blk.str = "Run",
  trt.str = "ANI",
  decimal = TRUE,
  digits = 4
)

summaryAovOnePhase(
  design.df = des,
  blk.str = "Run",
  trt.str = "Cag"
)

summaryAovTwoPhase(
  design.df = des,
  blk.str2 = "Run",
  blk.str1 = "ANI",
  trt.str = "Tag + Trt",
  decimal = TRUE,
  digits = 4
)$A

summaryAovTwoPhase(
  design.df = des,
  blk.str2 = "Run",
  blk.str1 = "Cag/ANI",
  trt.str = "Tag + Trt",
  decimal = TRUE,
  digits = 4
)$A

des$y = 1

str(des)

des$Run <- factor(des$Run)
des$Cag <- factor(des$Cag)
des$Tag <- factor(des$Tag)

summary(aov(y ~ Tag + Trt + Error(Run + Cag/ANI), des))

summary(aov(y ~ Cag/ANI + Error(Run), des))


for(i in list.files("OptimalDesigns/RCBD//")){
  des <- read.csv(paste0("OptimalDesigns/RCBD/", i))
  if(any(is.na(des$Tag)))
    print(i)
  
  
}
  library(optimTE)


des <- read.csv("OptimalDesigns/BIBD/designBIBD74724.csv")


out <- matrix(paste0(des$Cag), ncol = nlevels(factor(des$Tag)), 
              nrow = max(des$Run), byrow = TRUE)


summaryAovOnePhase(
  design.df = des,
  blk.str = "Cag/Ani",
  trt.str = "Trt",
  decimal = TRUE,
  digits = 4
)$A


summaryAovTwoPhase(
  design.df = des,
  blk.str2 = "Run",
  blk.str1 = "Cag/Ani",
  trt.str = "Tag + Trt",
  decimal = TRUE,
  digits = 4
)


des <- read.csv("OptimalDesigns/BIBD/designBIBD65624.csv")


out <- matrix(paste0(des$Cag), ncol = nlevels(factor(des$Tag)), 
              nrow = max(des$Run), byrow = TRUE)
out

summaryAovOnePhase(
  design.df = des,
  blk.str = "Cag/Ani",
  trt.str = "Trt",
  decimal = TRUE,
  digits = 4
)$A

summaryAovTwoPhase(
  design.df = des,
  blk.str2 = "Run",
  blk.str1 = "Cag/Ani",
  trt.str = "Tag + Trt",
  decimal = TRUE,
  digits = 4
)



















