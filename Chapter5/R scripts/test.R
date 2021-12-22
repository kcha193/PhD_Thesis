
des1 <- 
data.frame(Block = rep(1:4, each = 3), 
           Unit = LETTERS[c(2,1,1,3,2,2,4,3,3,1,4,4)],
           Trt = letters[c(2,1,1,1,2,2,2,1,1,1,2,2)])

pB <- table(1:12, des1$Block)
pU <- table(1:12, des1$Unit)
pT <- table(1:12, des1$Trt)

eigen(projMat(pU) %*% (projMat(pB) - K(12)) %*% projMat(pU))

fractions(eigen(projMat(pU) %*% (identityMat(12)- projMat(pB)) %*% projMat(pU))$val[1:3])

eigen(projMat(pT) %*% projMat(pU) %*% (projMat(pB) - K(12)) %*% 
        projMat(pU) %*% projMat(pT))$vec[,1]

round(eigen(projMat(pT) %*% projMat(pU) %*% (identityMat(12)- projMat(pB))  
      %*% projMat(pU) %*% projMat(pT))$vec[,1], 3)


eigen(projMat(pT)  %*% (projMat(pU) - K(12)) %*%  projMat(pT))

eigen(projMat(pT) %*% projMat(pU) %*% (identityMat(12)- projMat(pB))  
      %*% projMat(pU) %*% projMat(pT))


summaryAovOnePhase(des1, trt.str = "Unit", blk.str = "Block")

summaryAovTwoPhase(des1, trt.str = "Trt", blk.str1 = "Unit", 
                   blk.str2 = "Block", decimal = TRUE, latex = TRUE)



des1 <- 
  data.frame(Block = rep(1:4, each = 3), 
             Unit =  LETTERS[c(1,2,3,4,1,2,3,4,1,2,3,4)],
             Trt = c(1,2,1,2,1,2,1,2,1,2,1,2))

summaryAovOnePhase(des1, trt.str = "Unit", blk.str = "Block")

summaryAovTwoPhase(des1, trt.str = "Trt", blk.str1 = "Unit", blk.str2 = "Block")





des1 <- 
  data.frame(Block = rep(1:4, each = 3), 
             Unit =  LETTERS[c(1,2,3,4,1,2,3,4,1,2,3,4)],
             Trt = c(1,2,3,1,2,3,1,2,3,1,2,3))

summaryAovOnePhase(des1, trt.str = "Unit", blk.str = "Block")

summaryAovTwoPhase(des1, trt.str = "Trt", blk.str1 = "Unit", blk.str2 = "Block")








designPhase2 <- local({
  
  #Phase 2 Block structure
  nOc = 2; nJu = 6; nIn = 3; nSt = 4; nPos = 4
  Occ  <- factor(rep(1:nOc, each = nJu*nIn*nSt*nPos))
  Jud  <- factor(rep(1:nJu, times = nOc, each = nIn*nSt*nPos))
  Int  <- factor(rep(1:nIn, times = nOc*nJu,  each = nSt*nPos))
  Sit  <- factor(rep(1:nSt, times = nOc*nJu*nIn,  each = nPos))
  Pos <- factor(rep(1:nPos, times = nOc*nJu*nIn*nSt))
  
  #Phase 1 Block structure
  nR = 3; nS = 2; nC = 4; nH = 2
  RowSeq = c(1,3,2, 2,1,3, 3,2,1)
  Row      <- factor(rep(RowSeq, times = nC, each = nSt*nPos))
  
  Squ      <- Occ
  
  ColSeq1 = c(1,3,2,4, 1,2,3,4, 1,2,4,3)
  ColSeq2 = c(2,4,1,3, 3,4,1,2, 4,3,1,2)
  
  Col      <- factor(rep(c(rep(ColSeq1, times = nIn, each = nPos),
                           rep(ColSeq2, times = nIn, each = nPos)), times = nOc))
  Hal      <- factor(rep(1:nH, times = nOc*nJu*nIn*nSt, each = nH ))
  
  #Phase 1 Treatment structure
  Tre = factor(rep(c( 4, 2, 1, 3,  2, 3, 4, 1,  1, 2, 4, 3,
                      1, 3, 2, 4,  4, 1, 2, 3,  2, 3, 1, 4,
                      2, 4, 3, 1,  1, 2, 3, 4,  4, 1, 3, 2,
                      1, 3, 4, 2,  4, 1, 2, 3,  4, 3, 1, 2,
                      2, 4, 1, 3,  2, 3, 4, 1,  1, 4, 2, 3,
                      3, 1, 2, 4,  3, 4, 1, 2,  3, 2, 4, 1,
                      
                      2, 4, 1, 3,  4, 3, 2, 1,  3, 2, 4, 1,
                      3, 1, 2, 4,  2, 1, 4, 3,  4, 3, 1, 2,
                      4, 2, 3, 1,  3, 2, 1, 4,  2, 1, 3, 4,
                      1, 3, 2, 4,  2, 1, 4, 3,  4, 1, 3, 2,
                      2, 4, 3, 1,  4, 3, 2, 1,  1, 2, 4, 3,
                      3, 1, 4, 2,  1, 4, 3, 2,  3, 4, 2, 1), each = nPos))
  
  Met <- factor(LETTERS[Hal])
  
  data.frame(Tre, Met,  Hal, Col, Squ, Row, Pos, Jud, Sit, Int, Occ)
})
designPhase2

summaryAovTwoPhase(designPhase2,   blk.str2 = "(Occ/Int/Sit)*Jud", 
                   blk.str1 = "(Row*(Squ/Col))/Hal",
                   trt.str = "Tre")



summaryAovTwoPhase(designPhase2,  blk.str1 = "(Row*(Squ/Col))/Hal",
                   blk.str2 = "((Occ/Int/Sit)*Jud)/Pos", 
                   trt.str = "Tre*Met", decimal = TRUE, 
                   digits = 2, table.legend = TRUE)

