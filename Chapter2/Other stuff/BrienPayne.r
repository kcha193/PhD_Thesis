#BrienPayne Example
#Phase 1 Field phase

library(infoDecompuTE)


designPhase1 <- local({

   nR = 3; nS = 2; nC = 4; nH = 2

   Row      <- factor(rep(1:nR, each = nS*nC*nH))
   Squ      <- factor(rep(1:nS, times = nR, each = nC*nH))
   Col      <- factor(rep(1:nC, times = nR*nS, each = nH))
   Hal      <- factor(rep(1:nH, times = nR*nS*nC))

   Row.Squ  <- interaction(Row, Squ)
   Squ.Col  <- interaction(Squ, Col)
   Row.Squ.Col <- interaction(Row, Squ, Col)
   Row.Squ.Col.Hal <- interaction(Row, Squ, Col, Hal)

   Tre <-  factor(rep(c( 4, 1, 2, 3, 2, 1, 4, 3,
                         1, 2, 3, 4, 3, 2, 1, 4,
                         2, 3, 4, 1, 4, 3, 2, 1), each = 2))
   Pos <- Hal
   Inter <- interaction(Tre, Pos)

   #data.frame(Tre, Pos, Inter, Row.Squ.Col.Hal, Row.Squ.Col, Squ.Col, Row.Squ, Squ, Row)
   data.frame(Tre, Pos, Inter, Row, Squ, Col, Hal)
})
designPhase1


summaryAovOnePhase(designPhase1, blk.str = "(Row*(Squ/Col))/Hal", trt.str = "Tre",
table.legend = TRUE)



summaryAovOnePhase(designPhase1, blk.str = "(Row*(Squ/Col))/Hal", trt.str = "Tre")

#BrienPayne Example
#Phase 2 wine tasting phase

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

old <- proc.time()
summaryAovTwoPhase(designPhase2,  blk.str1 = "(Row*(Squ/Col))/Hal",
                                  blk.str2 = "((Occ/Int/Sit)*Jud)/Pos", 
                                  trt.str = "Tre*Met", decimal = TRUE, 
								  digits = 2, table.legend = TRUE)

proc.time() - old



###############################################################################


design = read.csv(file.choose(), header = TRUE, colClasses = c(rep("factor", 11), "numeric"))


summaryAovTwoPhase(design,  blk.str1 = "(Rows*(Squares/Columns))/Halfplot",
                                  blk.str2 = "((Occasion/Interval/Sittings)*Judges)/Position", 
                                  trt.str = "Trellis*Method", table.legend = TRUE)
  
 summaryAovOnePhase(design,  trt.str = "(Rows*(Squares/Columns))/Halfplot",
                                  blk.str = "((Occasion/Interval/Sittings)*Judges)/Position", table.legend = TRUE)
                                
summaryAovTwoPhase(design,  blk.str1 = "(Rows*(Squares/Columns))/Halfplot",
                                  blk.str2 = "((Occasion/Interval/Sittings)*Judges)/Position", 
                                  trt.str = "Trellis*Method",
                                  response = design$Score, table.legend = TRUE)


summaryAovTwoPhase(design,  blk.str1 = "(Rows*(Squares/Columns))/Halfplot",
                                  blk.str2 = "((Occasion/Interval/Sittings)*Judges)/Position", 
                                  trt.str = "Trellis*Method",
                                  response = design$Score, latex = TRUE)

