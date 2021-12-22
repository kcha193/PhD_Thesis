

summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Trt + Tag")

summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")

blk.contr = test(X.trt = Ani.mat,   (mI(n) - Pb) %*%  (mI(n) - Pb1), Rep=n/ncol(Ani.mat))$e.vec

blk.contr = test(X.trt = Ani.mat,   (mI(n) - Pb), Rep=n/ncol(Ani.mat))$e.vec

blk.contr =eigen(t(Ani.mat) %*%   (mI(n) - Pb) %*% Ani.mat)$vec

#blk.contr =svd(t(Ani.mat) %*%   (mI(n) - Pb) %*% Ani.mat)$v

A1 = blk.contr[,1][design.df$Ani]
A2 = blk.contr[,2][design.df$Ani]
A3 = blk.contr[,3][design.df$Ani]
A4 = blk.contr[,4][design.df$Ani]
A5 = blk.contr[,5][design.df$Ani]
A6 = blk.contr[,6][design.df$Ani]
A7 = blk.contr[,7][design.df$Ani]
A8 = blk.contr[,8][design.df$Ani]
A9 = blk.contr[,9][design.df$Ani]
A10 = blk.contr[,10][design.df$Ani]
A11 = blk.contr[,11][design.df$Ani]
A12 = blk.contr[,12][design.df$Ani]
A13 = blk.contr[,13][design.df$Ani]
A14 = blk.contr[,14][design.df$Ani]
A15 = blk.contr[,15][design.df$Ani]

blk.contr= list(Ani = list(A1 = A1, A2 = A2, A3 = A3, A4 = A4, A5 = A5, A6 = A6,
A7 = A7, A8 = A8, A9= A9, A10 = A10, A11 = A11, A12 = A12, A13 = A13, A14 = A14, A15 = A15))


summary.aov.onePhase(design.df,  blk.str = "Run", trt.str = "Ani", trt.contr = blk.contr)


blk.contr = test(X.trt = Ani.mat,  blk.proj, Rep=n/ncol(Ani.mat))$e.vec

blk.contr= list(Ani = list(A1 = cbind(A1, A2, A3), A2 = A4, A3 = cbind(A5, A6), A4 = cbind(A7, A8, A9), A5 = A10))


blk.contr= list(Ani = list(A1 = cbind(A1, A2, A3, A4, A5), 
                        A2 = cbind(A6, A7, A8, A9)))


summary.aov.onePhase(design.df,  blk.str = "Run", trt.str = "Ani", trt.contr = blk.contr)

summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", 
      trt.str = "Tag + Trt", blk.contr = blk.contr)



A1 = c(1,1,1,-1, -1, -1)
A2 = c(1,-1, 0, 1, -1, 0)
A3 = c(-1, -1 , 2, -1, -1, 2)
A4 = A1 * A2
A5 = A1 * A3

A1 = A1[design.df$Ani]
A2 = A2[design.df$Ani]
A3 = A3[design.df$Ani]
A4 = A4[design.df$Ani]
A5 = A5[design.df$Ani]



