bestInd = scan()
1
9
6
10
8
4
3
7
2
11
5
12

new.Z1.mat =  Z1.mat[bestInd,]
  #new.Z1.mat=  Z1.mat

################################################################################
  #Set the design from the search

  colnames(new.Z1.mat) = sort(levels(interaction(newLETTERS[1:nZ1])))
  new.Z1.des = apply(new.Z1.mat, 1,  function(x)
  colnames(new.Z1.mat)[which(as.logical(x))])
  new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))

  design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

  design.df = cbind(design.df,
              Ani = t(new.Z1.des),
              Trt = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                            as.character(phase1DesignEX1$Ani)),]$Trt)

################################################################################
  #Theortical ANOVA table


summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")

t(new.Z1.mat) %*% Zb
t(new.Z1.mat) %*% Zt


test(X.trt = Z1.mat[bestInd,], Pb - mK(n) , Rep=Z1.rep)
test(X.trt = Z1.mat[bestInd,], Pb1 - mK(n) , Rep=Z1.rep)

info.mat =  matMulti(Pb - mK(n) , Z1.mat[bestInd,], C.ani)


C.ani = mI(nAni) - mK(nAni)

info.mat =  matMulti(Pb - mK(n) , Z1.mat[bestInd,], C.ani)
blk.contr = test(X.trt = mI(nAni), info.mat , Rep=Z1.rep)$e.vec

blk.contr1 = blk.contr[,1][as.numeric(design.df$Ani)]
blk.contr2 = blk.contr[,2][as.numeric(design.df$Ani)]


info.mat = matMulti(Pb1 - mK(n) , Z1.mat[bestInd,], C.ani)
info.mat
fractions(info.mat)
test(X.trt = mI(nAni), info.mat , Rep=Z1.rep)
blk.contr = test(X.trt = mI(nAni), info.mat , Rep=Z1.rep)$e.vec
fractions(blk.contr)

blk.contr3 = blk.contr[,1][as.numeric(design.df$Ani)]

info.mat =  matMulti( blk.proj , Z1.mat[bestInd,], C.ani)
info.mat
fractions(info.mat)

test(X.trt = mI(nAni), info.mat , Rep=Z1.rep)

blk.contr = test(X.trt = mI(nAni), info.mat , Rep=Z1.rep)$e.vec
fractions(blk.contr)

blk.contr4 = blk.contr[,1][as.numeric(design.df$Ani)]
blk.contr5 = blk.contr[,2][as.numeric(design.df$Ani)]



summary.aov.onePhase(design.df, blk.str = "Run",  trt.str =  "Ani",
                      trt.contr = list(Ani = list(  blk1 = blk.contr1,
                                                    blk2 = blk.contr2,
                                                    blk3 = blk.contr3,
                                                    blk4 = blk.contr4,
                                                    blk5 = blk.contr5)))


summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt",
                      blk.contr = list(Ani = list(  blk1 = blk.contr1,
                                                    blk2 = blk.contr2,
                                                    blk3 = blk.contr3,
                                                    blk4 = blk.contr4,
                                                    blk5 = blk.contr5)))

summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag",
                      blk.contr = list(Ani = list(  blk1 = blk.contr1,
                                                    blk2 = blk.contr2,
                                                    blk3 = blk.contr3,
                                                    blk4 = blk.contr4,
                                                    blk5 = blk.contr5)))



info.mat =  matMulti((mI(n) - Pb)  , Z1.mat[bestInd,], C.ani)

test(X.trt = Zt, matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]), Rep= n/4 )


t(Zt) %*% matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]) %*% Zt

fractions((t(Zt) %*% matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]) %*% Zt) / (n/4))


tag.contr = test(X.trt = Zt, matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]), Rep= n/4 )$e.vec

test(X.trt = Zt, matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]), Rep= n/4 )

tag.contr1 = tag.contr[,1][as.numeric(design.df$Tag)]



t(Zt) %*% ((mI(n) - Pb) - matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,])) %*% Zt

fractions((t(Zt) %*% ((mI(n) - Pb) - matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,])) %*% Zt) / (n/4))

test(X.trt = Zt,  (mI(n) - Pb) - matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]), Rep= n/4 )

tag.contr = test(X.trt = Zt,  (mI(n) - Pb) - matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]), Rep= n/4 )$e.vec
fractions(tag.contr)

tag.contr2 = tag.contr[,1][as.numeric(design.df$Tag)]
tag.contr3 = tag.contr[,2][as.numeric(design.df$Tag)]



summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt",
                      trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3),
                                       Trt = NA))



summary.aov.onePhase(design.df, blk.str = "Run",  trt.str = "Tag + Trt",
                      trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3),
                                       Trt = NA))

summary.aov.onePhase(design.df, blk.str = "Ani",  trt.str = "Trt")


summary.aov.onePhase(design.df, blk.str = "Ani",  trt.str = "Tag + Trt",
                      trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3),
                                       Trt = NA))


summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt",
                      trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3),
                                       Trt = NA))

summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt",
                      blk.contr = list(Ani = list(  blk1 = blk.contr1,
                                                    blk2 = blk.contr2,
                                                    blk3 = blk.contr3,
                                                    blk4 = blk.contr4,
                                                    blk5 = blk.contr5)))

summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag",
                      blk.contr = list(Ani = list(  blk1 = blk.contr1,
                                                    blk2 = blk.contr2,
                                                    blk3 = blk.contr3,
                                                    blk4 = blk.contr4,
                                                    blk5 = blk.contr5)))

summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag",
                      blk.contr = list(Ani = list(blk1 = blk.contr1,
                                        blk2 = blk.contr2,
                                        blk3 = blk.contr3,
                                        blk4 = blk.contr4,
                                        blk5 = blk.contr5)),
                      trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3)))

summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt",
                      blk.contr = list(Ani = list(blk1 = blk.contr1,
                                        blk2 = blk.contr2,
                                        blk3 = blk.contr3,
                                        blk4 = blk.contr4,
                                        blk5 = blk.contr5)),
                      trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3),
                                      Trt = NA))

summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")


Ani.contr = cbind(blk.contr1, blk.contr2, blk.contr3, blk.contr4, blk.contr5)
Tag.contr = cbind(tag.contr1, tag.contr2, tag.contr3)

t(Ani.contr )  %*% Tag.contr

round(t(Ani.contr )  %*% Tag.contr, 5)
