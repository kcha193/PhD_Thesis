#21/05/2013 9:02:15 a.m.
#Computing the stratum variance using the REML based method from Payne and Tobias  (1992)



WiBlk = function(design.df) {
    
    n = nrow(design.df)    
       
    nBlk = nlevels(design.df$Run)
    nPlot = nlevels(design.df$Tag)
    
    nAni = nlevels(design.df$Ani)
    Run.mat = with(design.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
    Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
    Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))
    
    Trt.mat = Trt.mat[, sort(colnames(Trt.mat))]

     C.ani = (identityMat(nAni) - K(nAni))
     
    Pb = projMat(Run.mat)
    
    blk.proj = (identityMat(n) - Pb) 
    
    
    PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani, Ani.mat)
    
    return(PP1)
}

BwBlk = function(design.df) {
    
    n = nrow(design.df)    
       
    nBlk = nlevels(design.df$Run)
    nPlot = nlevels(design.df$Tag)
    
    nAni = nlevels(design.df$Ani)
    Run.mat = with(design.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
    Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
    Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))
    
    Trt.mat = Trt.mat[, sort(colnames(Trt.mat))]

     C.ani = (identityMat(nAni) - K(nAni))
     
    Pb = projMat(Run.mat)
    
    blk.proj = ( Pb - K(n)) 
    
    
    PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani, Ani.mat)
    
    return(PP1)
}



trtProj = function(design.df) {
    
   n = nrow(design.df)    
       
    nBlk = nlevels(design.df$Run)
    nPlot = nlevels(design.df$Tag)
    
    nAni = nlevels(design.df$Ani)
    Run.mat = with(design.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
    Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
    Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))
    
    Trt.mat = Trt.mat[, sort(colnames(Trt.mat))]

     C.ani = (identityMat(nAni) - K(nAni))
     
    Pb = projMat(Run.mat)
    
    blk.proj = (identityMat(n) - Pb) 
    
    
    PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani, Ani.mat)

    return(ginv(t(Trt.mat) %*% PP1 %*% Trt.mat) %*% t(Trt.mat) %*% PP1)
}


summaryAovTwoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")

gamma.run = 0.001
gamma.ani = 1
VC.resid = 1

true.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.run * VC.resid))

run.eff = rnorm(nlevels(design$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
ani.eff = rnorm(nlevels(design$Ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
trt.eff = runif(nlevels(design$Trt), 0, 2)
tag.eff = runif(nlevels(design$Tag), 0, 1)
res.eff = rnorm(nrow(design), mean = 0, sd = sqrt(VC.resid))
                
y = gm + with(design, run.eff[Run] + ani.eff[Ani] +  tag.eff[Tag] + trt.eff[Trt]) + res.eff

summaryAovTwoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", 
trt.str = "Tag + Trt", response = y)

Sa1 = BwBlk(design.df)
Sa2 = WiBlk(design.df)

trtEff = trtProj(design.df) %*% y

Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))
 
X = cbind(1, Trt.mat, Tag.mat)
L = as.numeric( X %*% ginv(t(X) %*% X) %*% t(X) %*% y)

Sa = projMat(Ani.mat)

MLEstraVar1 = t(y - L) %*% Sa1 %*% (y - L)/tr(Sa1)    # downward estimate
MLEstraVar2 = t(y - L) %*% Sa2 %*% (y - L)/tr(Sa2)    # downward estimate
     
w2 = ((8/9)/MLEstraVar2)/(MLEstraVar1 + (8/9)/MLEstraVar2)
w1 = MLEstraVar1/(MLEstraVar1 + (8/9)/MLEstraVar2)

d2 = tr(Sa2) -  w2 *tr(projMat(Trt.mat ))
d1 = tr(Sa1) -  w1 *tr(projMat(Trt.mat ))

MLEstraVar1 = t(y - L) %*% Sa1 %*% (y - L)/d1    # downward estimate
MLEstraVar2 = t(y - L) %*% Sa2 %*% (y - L)/d2   # downward estimate

c(MLEstraVar1, MLEstraVar2)

(d2/MLEstraVar2 + d1/MLEstraVar1)^2/((1/MLEstraVar1)^2/d1 + (1/MLEstraVar2)^2/d2)
      
  (MLEstraVar2 )^2/(MLEstraVar1^2/d1 + MLEstraVar2^2/d2)
  
(MLEstraVar2/(6 * 8/9) + MLEstraVar1/6)^2 /((MLEstraVar1/6)^2/1 + (MLEstraVar2/(6 * 8/9))^2/2   )
  
      
aov.table = summaryAovTwoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", 
 trt.str = "Tag + Trt", response = y)
                    
getVcEDF(aov.table = aov.table, row.MS = NA, true.VC =  true.VC, neg.VC = FALSE)


####################################################################################################


WiBlk = function(design1.df) {
    
    n = nrow(design1.df)   
    
    
    nBlk = nlevels(design1.df$Run)
    nPlot = nlevels(design1.df$Tag)
    nCag = nlevels(design1.df$Cag)
    
    nAni = nlevels(design1.df$Ani)
    Run.mat = with(design1.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design1.df, as.matrix(table(1:n, Tag)))
    Cag.mat = with(design1.df, as.matrix(table(1:n, Cag)))
    Ani.mat = with(design1.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
        Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))
    
    Trt.mat = Trt.mat[, sort(colnames(Trt.mat))]

   
    C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
    C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
    cage.Rep = n/nCag
    
    C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
    
    Pb = projMat(Run.mat)
    Pb1 = projMat(Tag.mat)
    
    blk.proj = (identityMat(n) - Pb)
    
    info.mat = matMulti(blk.proj, Ani.mat, C.cage)
    
    if (any(abs(info.mat) > 1e-07)) {
        PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)), C.cage, Ani.mat)
        
        PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
    } else {
        PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani, Ani.mat)
        
    }
    return(PP1)
}


BwBlk = function(design1.df) {
    
    n = nrow(design1.df)   
    
    
    nBlk = nlevels(design1.df$Run)
    nPlot = nlevels(design1.df$Tag)
    nCag = nlevels(design1.df$Cag)
    
    nAni = nlevels(design1.df$Ani)
    Run.mat = with(design1.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design1.df, as.matrix(table(1:n, Tag)))
    Cag.mat = with(design1.df, as.matrix(table(1:n, Cag)))
    Ani.mat = with(design1.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
        Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))
    
    Trt.mat = Trt.mat[, sort(colnames(Trt.mat))]

   
    C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
    C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
    cage.Rep = n/nCag
    
    C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
    
    Pb = projMat(Run.mat)
    Pb1 = projMat(Tag.mat)
    
    blk.proj = (Pb - K(n))
    
    info.mat = matMulti(blk.proj, Ani.mat, C.cage)
    
    if (any(abs(info.mat) > 1e-07)) {
        PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)), C.cage, Ani.mat)
        
        PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
    } else {
        PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani, Ani.mat)
        
    }
    return(PP1)
}


basicContr = 
function(design.df){

    n = nrow(design.df)
    nBlk = nlevels(design.df$Run)
    nPlot = nlevels(design.df$Tag)
    nCag = nlevels(design.df$Cag)

    nAni = nlevels(design.df$Ani)
    Run.mat = with(design.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
    Cag.mat = with(design.df, as.matrix(table(1:n, Cag)))
    Ani.mat = with(design.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
    Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))
    
    Trt.mat = Trt.mat[, sort(colnames(Trt.mat))]

    C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
    C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
    cage.Rep = n/nCag

    C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))

    Pb = projMat(Run.mat)
    Pb1 = projMat(Tag.mat)

    blk.proj = (identityMat(n) - Pb) %*%  (identityMat(n) - Pb1)

     
      info.mat = matMulti(blk.proj, Ani.mat, C.cage)
  
      if (any(abs(info.mat) > 1e-07)) {
          PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)),
              C.cage, Ani.mat)
  
          PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
      } else {
          PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani,
              Ani.mat)
  
      }
  
      test.RBD(X.trt = Trt.mat, PP1, Rep = n/ncol(Trt.mat))$e.vec

}


summaryAovTwoPhase(design6, blk.str2 = "Run", blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")



trt.contr = basicContr(design6)

trt.contr1 = trt.contr[,1][as.numeric(design6$Trt)]
trt.contr2 = trt.contr[,2][as.numeric(design6$Trt)]
trt.contr3 = trt.contr[,3][as.numeric(design6$Trt)]
trt.contr4 = trt.contr[,4][as.numeric(design6$Trt)]
trt.contr5 = trt.contr[,5][as.numeric(design6$Trt)]


summaryAovTwoPhase(design6,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt",trt.contr = list(Tag = NA,
                                       Trt = list(  "1" = trt.contr1,
                                                    "2" = trt.contr2,
                                                    "3" = trt.contr3,
                                                    "4" = trt.contr4,
                                                    "5" = trt.contr5)))


summaryAovTwoPhase(design6,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt",trt.contr = list(Tag = NA,
                                       Trt = list("1" = trt.contr1,
                                                  "3" = trt.contr3)))


summaryAovTwoPhase(design6,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt",trt.contr = list(Tag = NA,
                                       Trt = list( "2" = trt.contr2)))

summaryAovTwoPhase(design6,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt",trt.contr = list(Tag = NA,
                                       Trt = list( "4" = trt.contr4)))

summaryAovTwoPhase(design6,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt",trt.contr = list(Tag = NA,
                                       Trt = list( "5" = trt.contr5)))


gamma.run = 0 
gamma.cag = 10
gamma.ani = 10       
VC.resid = 1        

true.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag * VC.resid), (gamma.run * VC.resid))

run.eff = rnorm(nlevels(design6$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
cag.eff = rnorm(nlevels(design6$Cag), mean = 0, sd = sqrt(gamma.cag * VC.resid))
ani.eff = rnorm(nlevels(design6$ani), mean = 0, sd = sqrt(gamma.ani* VC.resid))
trt.eff = runif(nlevels(design6$Trt), 0, 10)
tag.eff = runif(nlevels(design6$Tag), 0, 1)
res.eff = rnorm(nrow(design6), mean = 0, sd = sqrt(VC.resid))
gm = 10

y6 = gm + with(design6, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + tag.eff[Tag] + 
    trt.eff[Trt]) + res.eff

aov.table = summaryAovTwoPhase(design6,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt",trt.contr = list(Tag = NA,
                                       Trt = list( "2" = trt.contr2)), response = y6)

getVcEDF(aov.table = aov.table, row.MS = NA, true.VC =  true.VC, neg.VC = FALSE)

Sa1 = BwBlk(design6)
Sa2 = WiBlk(design6)


Tag.mat = with(design6, as.matrix(table(1:n, Tag)))
Trt.mat = with(design6, as.matrix(table(1:n, Trt)))
 
X = cbind(1, t(trt.contr[,2] %*% t(trt.contr[,2]) %*% t(Trt.mat)), Tag.mat)

L = as.numeric( projMat(X) %*% y6)

MLEstraVar1 = t(y6 - L) %*% Sa1 %*% (y6 - L)/(tr(Sa1))      # downward estimate
MLEstraVar2 = t(y6 - L) %*% Sa2 %*% (y6 - L)/(tr(Sa2))    # downward estimate
     
w2 = ((6 * 297200/332313)/MLEstraVar2)/(((6 * 75658/716035)/MLEstraVar1) + ((6 * 297200/332313)/MLEstraVar2))
w1 = ((6 * 75658/716035)/MLEstraVar1)/(((6 * 75658/716035)/MLEstraVar1) + ((6 * 297200/332313)/MLEstraVar2))

d2 = tr(Sa2) -  w2 * 2
d1 = tr(Sa1) -  w1 *1 

MLEstraVar1 = t(y6 - L) %*% Sa1 %*% (y6 - L)/(tr(Sa1) - d1)      # downward estimate
MLEstraVar2 = t(y6 - L) %*% Sa2 %*% (y6 - L)/(tr(Sa2) - d2)    # downward estimate
     

c(MLEstraVar1, MLEstraVar2)



((MLEstraVar2/(8/9)) +  MLEstraVar1)^2/((MLEstraVar1)^2/d1 + (MLEstraVar2/(8/9))^2/d2)
   










