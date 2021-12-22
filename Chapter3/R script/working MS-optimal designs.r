bestInd = scan()
4
10
8
11
5
1
12
9
6
3
7
2



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


test(X.trt = Z1.mat[bestInd,], Pb - mK(n) , Rep=Z1.rep)
test(X.trt = Z1.mat[bestInd,], mI(n) - Pb , Rep=Z1.rep)

info.mat =  matMulti(mI(n) - Pb, Z1.mat[bestInd,], C.ani)


C.ani = mI(nAni) - mK(nAni)

info.mat =  matMulti(blk.proj , Z1.mat[bestInd,], C.ani)


blk.contr = test(X.trt = Z1.mat[bestInd,], blk.proj, Rep=Z1.rep)$e.vec

blk.contr = test(X.trt = mI(nAni), info.mat, Rep =Z1.rep)$e.vec

fractions(blk.contr)

blk.contr1 = blk.contr[,1][as.numeric(design.df$Ani)]
blk.contr2 = blk.contr[,2][as.numeric(design.df$Ani)]
blk.contr3 = blk.contr[,3][as.numeric(design.df$Ani)]
blk.contr4 = blk.contr[,4][as.numeric(design.df$Ani)]
blk.contr5 = blk.contr[,5][as.numeric(design.df$Ani)]


summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt",
                      blk.contr = list(Ani = list(  blk1 = blk.contr1,
                                                    blk2 = blk.contr2,
                                                    blk3 = blk.contr3,
                                                    blk4 = blk.contr4,
                                                    blk5 = blk.contr5)))


summary.aov.onePhase(design.df, blk.str = "Run",  trt.str =  "Ani",
                      trt.contr = list(Ani = list(  blk1 = blk.contr1,
                                                    blk2 = blk.contr2,
                                                    blk3 = blk.contr3,
                                                    blk4 = blk.contr4,
                                                    blk5 = blk.contr5)))


info.mat =  matMulti((mI(n) - Pb)  , Z1.mat[bestInd,], C.ani)

test(X.trt = Zt, matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]), Rep= n/4 )


t(Zt) %*% matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]) %*% Zt

fractions(t(Zt) %*% matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]) %*% Zt)


tag.contr = test(X.trt = Zt, matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]), Rep= n/4 )$e.vec

tag.contr1 = tag.contr[,1][as.numeric(design.df$Tag)]
tag.contr2 = tag.contr[,2][as.numeric(design.df$Tag)]
tag.contr3 = tag.contr[,3][as.numeric(design.df$Tag)]

summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt",
                      trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3),
                                       Trt = NA))

test(X.trt = Zt,  (mI(n) - Pb) - matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]), Rep= n/4 )

tag.contr = test(X.trt = Zt,  (mI(n) - Pb) - matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,]), Rep= n/4 )$e.vec


t(Zt) %*% ((mI(n) - Pb) - matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,])) %*% Zt

fractions(t(Zt) %*% ((mI(n) - Pb) - matMulti1((mI(n) - Pb) , ginv( info.mat), C.ani, Z1.mat[bestInd,])) %*% Zt)


tag.contr1 = tag.contr[,1][as.numeric(design.df$Tag)]
tag.contr2 = tag.contr[,2][as.numeric(design.df$Tag)]
tag.contr3 = tag.contr[,3][as.numeric(design.df$Tag)]


summary.aov.onePhase(design.df, blk.str = "Run",  trt.str = "Tag + Trt",
                      trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3),
                                       Trt = NA))

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

summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag",
                      blk.contr = list(Ani = list(  blk1 = blk.contr1,
                                                    blk2 = blk.contr2,
                                                    blk3 = blk.contr3,
                                                    blk4 = blk.contr4,
                                                    blk5 = blk.contr5)))

summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt",
                      blk.contr = list(Ani = list( blk5 = blk.contr5,
                                          blk1 = blk.contr1,
                                        blk2 = blk.contr2,
                                        blk3 = blk.contr3,
                                        blk4 = blk.contr4
                                        )),
                      trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3),
                                       Trt = NA))

summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag",
                      blk.contr = list(Ani = list( blk5 = blk.contr5,
                                          blk1 = blk.contr1,
                                        blk2 = blk.contr2,
                                        blk3 = blk.contr3,
                                        blk4 = blk.contr4
                                        )),
                      trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3)))


summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")


summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt",
   blk.contr = list(Ani = list(blk1 = cbind(blk.contr1,blk.contr3,blk.contr4),
                               blk2 = cbind(blk.contr2, blk.contr5))),
   trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3),
                                       Trt = NA))


summary.aov.onePhase(design.df, blk.str = "Run",  trt.str = "Tag",
                      trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3)))

summary.aov.onePhase(design.df, blk.str = "Run",  trt.str = "Ani",
                      trt.contr = list(Ani = list(  blk5 = blk.contr5,
                                                    blk1 = blk.contr1,
                                                    blk2 = blk.contr2,
                                                    blk3 = blk.contr3,
                                                    blk4 = blk.contr4)))

summary.aov.onePhase(design.df, blk.str = "Run",  trt.str = "Ani + Tag",
                      trt.contr = list( Ani = list(  blk5 = blk.contr5,
                                                    blk1 = blk.contr1,
                                                    blk2 = blk.contr2,
                                                    blk3 = blk.contr3,
                                                    blk4 = blk.contr4),
                                        Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3)))

summary.aov.onePhase(design.df, blk.str = "Run",  trt.str = "Tag + Ani",
                      trt.contr = list(Tag = list(  "1" = tag.contr1,
                                                    "2" = tag.contr2,
                                                    "3" = tag.contr3),
                                        Ani = list(  blk5 = blk.contr5,
                                                    blk1 = blk.contr1,
                                                    blk2 = blk.contr2,
                                                    blk3 = blk.contr3,
                                                    blk4 = blk.contr4)))

Ani.contr = cbind(blk.contr1, blk.contr2, blk.contr3, blk.contr4, blk.contr5)
Tag.contr = cbind(tag.contr1, tag.contr2, tag.contr3)

t(Ani.contr )  %*% Tag.contr


