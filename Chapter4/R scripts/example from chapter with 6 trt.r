
pairwiseContr = function(design1.df){

n = nrow(design1.df)
Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))

states = lapply(1:ncol(Trt.mat), function(x) t(cbind(x, (x+1):ncol(Trt.mat))))
states = states[-length(states)]

states= t(matrix(c(states, recursive=TRUE), nrow = 2))

pairMat = matrix(0, nrow = ncol(Trt.mat), ncol = 1)
 
for( i in 1:nrow(states)){

pairMat = cbind(pairMat, pairMat[,1])

pairMat[states[i,1],i+1] = 0.5
pairMat[states[i,2],i+1] = -0.5

}

pairMat = pairMat[,-1]


nBlk = nlevels(design1.df$Run)
nPlot = nlevels(design1.df$Tag)
nCag = nlevels(design1.df$Cag)

nAni = nlevels(design1.df$Ani)
Run.mat = with(design1.df, as.matrix(table(1:n, Run)))
Tag.mat = with(design1.df, as.matrix(table(1:n, Tag)))
Cag.mat = with(design1.df, as.matrix(table(1:n, Cag)))
Ani.mat = with(design1.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))

C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
cage.Rep = n/nCag

C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))

Pb = projMat(Run.mat)
Pb1 = projMat(Tag.mat)

blk.proj = (identityMat(n) - Pb) %*% (identityMat(n) - Pb1)

info.mat = matMulti(blk.proj, Ani.mat, C.cage)

  if (any(abs(info.mat) > 1e-07)) {
      PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)),
          C.cage, Ani.mat)

      PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
  } else {
      PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani,
          Ani.mat)

  }

contr.e = (apply(pairMat, 2, function(c)  t(c) %*% ginv(t(Trt.mat) %*% Trt.mat) %*% c)/apply(pairMat, 2, function(c)  t(c) %*% ginv(t(Trt.mat) %*% PP1 %*% Trt.mat) %*% c))

#browser()
      
#names.contr.e = apply(states, 1, function(x) paste(x, collapse = "-"))      

contrastMat = matrix(0, ncol= ncol(Trt.mat), nrow =  ncol(Trt.mat))

for( i in 1:length(contr.e)){   
contrastMat[states[i,2], states[i,1]] = contr.e[i]
}

return(contrastMat)
}
######################################################################################

pairwiseContr = function(design1.df){

  n = nrow(design1.df)
  Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))
  
  states = lapply(1:ncol(Trt.mat), function(x) t(cbind(x, (x+1):ncol(Trt.mat))))
  states = states[-length(states)]
  
  states= t(matrix(c(states, recursive=TRUE), nrow = 2))
  
  pairMat = matrix(0, nrow = ncol(Trt.mat), ncol = 1)
   
  for( i in 1:nrow(states)){
  
    pairMat = cbind(pairMat, pairMat[,1])
    
    pairMat[states[i,1],i+1] = 0.5
    pairMat[states[i,2],i+1] = -0.5
  
  }
  
  
  return(pairMat[,-1])
}

contrCompare = function(design1.df, contr){
  n = nrow(design1.df)
  Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))
  
  Trt.mat = Trt.mat[,sort(colnames(Trt.mat))]
  
  states = lapply(1:ncol(Trt.mat), function(x) t(cbind(x, (x+1):ncol(Trt.mat))))
  states = states[-length(states)]
  
  states= t(matrix(c(states, recursive=TRUE), nrow = 2))
  
   
  pairMat = contr
  
  
  nBlk = nlevels(design1.df$Run)
  nPlot = nlevels(design1.df$Tag)
  nCag = nlevels(design1.df$Cag)
  
  nAni = nlevels(design1.df$Ani)
  Run.mat = with(design1.df, as.matrix(table(1:n, Run)))
  Tag.mat = with(design1.df, as.matrix(table(1:n, Tag)))
  Cag.mat = with(design1.df, as.matrix(table(1:n, Cag)))
  Ani.mat = with(design1.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
  
  C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
  C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
  cage.Rep = n/nCag
  
  C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
  
  Pb = projMat(Run.mat)
  Pb1 = projMat(Tag.mat)
  
  blk.proj = (identityMat(n) - Pb) %*% (identityMat(n) - Pb1)
  
  info.mat = matMulti(blk.proj, Ani.mat, C.cage)
  
    if (any(abs(info.mat) > 1e-07)) {
        PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)),
            C.cage, Ani.mat)
  
        PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
    } else {
        PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani,
            Ani.mat)
  
    }
  
  contr.e = (apply(pairMat, 2, function(c)  t(c) %*% ginv(t(Trt.mat) %*% Trt.mat) %*% c)/apply(pairMat, 2, function(c)  t(c) %*% ginv(t(Trt.mat) %*% PP1 %*% Trt.mat) %*% c))
  
   names(contr.e) =   apply(states, 1, function(x) paste(colnames(Trt.mat)[x], sep = "", collapse = "-"))
  
  return(contr.e)
}


trtProj = function(design1.df){

n = nrow(design1.df)
Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))



nBlk = nlevels(design1.df$Run)
nPlot = nlevels(design1.df$Tag)
nCag = nlevels(design1.df$Cag)

nAni = nlevels(design1.df$Ani)
Run.mat = with(design1.df, as.matrix(table(1:n, Run)))
Tag.mat = with(design1.df, as.matrix(table(1:n, Tag)))
Cag.mat = with(design1.df, as.matrix(table(1:n, Cag)))
Ani.mat = with(design1.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))

C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
cage.Rep = n/nCag

C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))

Pb = projMat(Run.mat)
Pb1 = projMat(Tag.mat)

blk.proj = (identityMat(n) - Pb) %*% (identityMat(n) - Pb1)

info.mat = matMulti(blk.proj, Ani.mat, C.cage)

  if (any(abs(info.mat) > 1e-07)) {
      PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)),
          C.cage, Ani.mat)

      PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
  } else {
      PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani,
          Ani.mat)

  }

return( ginv(t(Trt.mat) %*% PP1 %*% Trt.mat) %*% t(Trt.mat) %*% PP1)
}

####################################################################################
  
  
sourceDir <- function(path, trace = TRUE, ...) {
  library(MASS)
   library(inline)
  library(compiler)
  library(formatR)
  library(Rcpp)
  library(RcppArmadillo)
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
  }
}

sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R")
sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/optimTE/R")




old = proc.time()
designCRD = optCRD(nTrt  = 6, bRep  = 3, tRep  = 2, nPlot = 4, iter  = 1000)
proc.time() - old

n = nrow(design.df)
nBlk = nlevels(design.df$Run)
nPlot = nlevels(design.df$Tag)
Run.mat = with(design.df, as.matrix(table(1:n, Run)))
Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))
Pb = projMat(Run.mat)
Pb1 = projMat(Tag.mat)

(ave.eff = test.CRD(X.trt = Trt.mat, (diag(n) - Pb) %*% (diag(n) - Pb1), Rep = n/ncol(Trt.mat))$ave.eff )

design1.df = optRBD( nTrt = 7, bRep  = 6, nCag  = 3, tRep  = 2, nPlot = 4, iter  = 3000, 
upperValue = ave.eff, confoundCag = FALSE)

design.summary.RBD(design1.df, FALSE)

pairwiseCompare(design1.df) 
    
design2.df = optRBD( nTrt = 7, bRep  = 6, nCag  = 3, tRep  = 2, nPlot = 4, iter  = 3000, 
upperValue = ave.eff, confoundCag = FALSE, resDF = 22)

design.summary.RBD(design2.df, FALSE)

pairwiseCompare(design2.df)

 summary(pairwiseCompare(design1.df))
 summary(pairwiseCompare(design2.df))

 sum(pairwiseCompare(design2.df) > pairwiseCompare(design1.df))
 which(pairwiseCompare(design2.df) > pairwiseCompare(design1.df))

 #####################################################################################

old = proc.time()
designCRD = optCRD(nTrt  = 6, bRep  = 3, tRep  = 2, nPlot = 4, iter  = 1000)
proc.time() - old

design.summary.CRD(design.df)

n = nrow(design.df)
nBlk = nlevels(design.df$Run)
nPlot = nlevels(design.df$Tag)
Run.mat = with(design.df, as.matrix(table(1:n, Run)))
Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))
Pb = projMat(Run.mat)
Pb1 = projMat(Tag.mat)

(ave.eff = test.CRD(X.trt = Trt.mat, (diag(n) - Pb) %*% (diag(n) - Pb1), Rep = n/ncol(Trt.mat))$ave.eff )

design1.df = optRBD( nTrt = 6, bRep  = 3, nCag  = 3, tRep  = 2, nPlot = 4, iter  = 3000, 
upperValue = ave.eff, confoundCag = FALSE)


design1new.df = optRBD( nTrt = 6, bRep  = 3, nCag  = 3, tRep  = 2, nPlot = 4, iter  = 3000, 
upperValue = ave.eff, confoundCag = FALSE)

design.summary.RBD(design1.df, FALSE)
design.summary.RBD(design1new.df)


(t1 = pairwiseCompare(design1.df)) 
(t2 = pairwiseCompare(design1new.df)) 

all((sort(t1[lower.tri(t1)])- sort(t2[lower.tri(t2)])) < 1e-7)
      
design2.df = optRBD( nTrt = 6, bRep  = 3, nCag  = 3, tRep  = 2, nPlot = 4, iter  = 3000, 
upperValue = ave.eff, confoundCag = FALSE, resDF = 5)

design2new.df = optRBD( nTrt = 6, bRep  = 3, nCag  = 3, tRep  = 2, nPlot = 4, iter  = 3000, 
upperValue = ave.eff, confoundCag = FALSE, resDF = 5)

design.summary.RBD(design2.df, FALSE)
design.summary.RBD(design2new.df, FALSE)

pairwiseContr(design6)


#####################################################################


obj = function(ind, design, contr){
 temp = contrCompare(design,  contr[ind,])
 
 #return(-max(temp))
 sd(temp)
}

swap = function(ind, design, contr) sample(length(ind),length(ind))



round(basicContr(design6)[,1:5], 4)

contrCompare(design6, basicContr(design6)[,1:5])

1/mean(1/contrCompare(design6, basicContr(design6)[,1:5]))


round(basicContr(design5), 4)

contrCompare(design5, basicContr(design5)[,1:5])

1/mean(1/contrCompare(design5, basicContr(design5)[,1:5]))

ans6 = optim(1:6, obj, swap, method = "SANN", design= design6, contr = basicContr(design6)[,1:5],
             control = list(maxit = 10000, temp = 0.1, trace = TRUE,
                            REPORT = 100))

contrCompare(design6,  basicContr(design6)[ans6$par,1:5])

ans5 = optim(1:6, obj, swap, method = "SANN", design= design5, contr = basicContr(design5)[,1:5],
             control = list(maxit = 10000, temp = .1, trace = TRUE,
                            REPORT = 100))

contrCompare(design5, basicContr(design5)[ans5$par,1:5])


contr1 = c(1, -1, 0,0,0,0)
contr2 = c( 0,0,1, -1,0,0)
contr3 = c( 0,0,0,0,1, -1)
contr4 = c(0.5, 0.5, 0.5,0.5,-1,-1)
contr5 = c(1,1, -1,-1,0,0)
  
contr = cbind(contr1, contr2, contr3, contr4, contr5)/2  
contr

contrCompare(design6, contr)

1/mean(1/contrCompare(design6, contr))

contrCompare(design5, contr)

1/mean(1/contrCompare(design5, contr))

ans6 = optim(1:6, obj, swap, method = "SANN", design= design6, contr = contr,
             control = list(maxit = 10000, temp = 0.1, trace = TRUE,
                            REPORT = 100))
contr[ans6$par,]
contrCompare(design6,  contr[ans6$par,])
1/mean(1/contrCompare(design6,  contr[ans6$par,]))


ans5 = optim(1:6, obj, swap, method = "SANN", design= design5, contr = contr,
             control = list(maxit = 10000, temp = .1, trace = TRUE,
                            REPORT = 100))
contr[ans5$par,]
contrCompare(design5, contr[ans5$par,])
contrCompare(design6, contr[ans5$par,])

1/mean(1/contrCompare(design5, contr[ans5$par,]))


contr1 = c(1, 1, 1,-1,-1,-1)
contr2 =c(0.5, 0.5, -1,0.5,0.5,-1)
contr3 = c( 1, -1,0, 1, -1, 0)
contr4 = contr1 * contr2
contr5 = contr1 * contr3

contr = cbind(contr1, contr2, contr3, contr4, contr5)/2  
contr

contrCompare(design6, contr)

1/mean(1/contrCompare(design6, contr))

contrCompare(design5, contr)

1/mean(1/contrCompare(design5, contr))

ans6 = optim(1:6, obj, swap, method = "SANN", design= design6, contr = contr,
             control = list(maxit = 10000, temp = 0.1, trace = TRUE,
                            REPORT = 100))
contr[ans6$par,]
contrCompare(design6,  contr[ans6$par,])
1/mean(1/contrCompare(design6,  contr[ans6$par,]))


ans5 = optim(1:6, obj, swap, method = "SANN", design= design5, contr = contr,
             control = list(maxit = 10000, temp = .1, trace = TRUE,
                            REPORT = 100))
contr[ans5$par,]
contrCompare(design5, contr[ans5$par,])
1/mean(1/contrCompare(design5, contr[ans5$par,]))


contr.e = contrCompare(design6, pairwiseContr(design6))

contrastMat = matrix(0, ncol= 6, nrow = 6)

for( i in 1:length(contr.e)){   
  contrastMat[states[i,2], states[i,1]] = contr.e[i]
}

contrastMat

#Pick the best 3 pairwise contrasts 

contr1 = c(0, 0, 0,1,0,-1)
contr2 = c( 1,0,-1, 0,0,0)
contr3 = c( 0,1,0,0,-1, 0)
contr4 = c(-1, 0,-1,1,0,1)
contr5 = c(0.5,-1, 0.5,0.5,-1,0.5)
  
contr = cbind(contr1, contr2, contr3, contr4, contr5)/2  
contr

contrCompare(design6, contr)


contr.e = contrCompare(design5, pairwiseContr(design5))

contrastMat = matrix(0, ncol= 6, nrow = 6)

for( i in 1:length(contr.e)){   
contrastMat[states[i,2], states[i,1]] = contr.e[i]
}

contrastMat

#Pick the best 3 pairwise contrasts 

contr1 = c(1, 0, 0,0,0,-1)
contr2 = c(0, 0,1, -1,0,0)
contr3 = c( 0,1,0,0,-1, 0)
contr4 = c(1, 0,-1,-1,0,1)
contr5 = c(0.5,-1, 0.5,0.5,-1,0.5)
  
contr = cbind(contr1, contr2, contr3, contr4, contr5)/2  
contr

contrCompare(design5, contr)


####################################################################################
design.summary.CRD(designCRD, FALSE)

design.summary.RBD(design6, FALSE)

design.summary.RBD(design5, FALSE)



#####################################################################################
#6 Treatments simulation   by normal distribution
design6 = design1.df
design5 = design2.df
design.df = design1.df

design.df$Ani =  with(design.df, interaction(Cag, Ani))
design6$ani =  with(design6, interaction(Cag, Ani))
design5$ani =  with(design5, interaction(Cag, Ani))

(anovaTable6 = summaryAovTwoPhase(design.df = design6 ,blk.str2 = "Run", 
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)

(anovaTable5 = summaryAovTwoPhase(design.df = design5 ,blk.str2 = "Run", 
              blk.str1 = "Cag/Ani", trt.str = "Tag + Trt")$A)

contr6.e = contrCompare(design6, pairwiseContr(design6))
contr5.e = contrCompare(design5, pairwiseContr(design5))
    
plot(contr6.e, contr5.e)
   
trtEff6 = trtProj(design6)
trtEff5 = trtProj(design5)


old = proc.time()
VC.resid = 1

gamma.run =  10
gamma.cag =  5   
gamma.ani =  2 #  c(10^((-4:4)/2))

nSim = 10000

temp6 = data.frame()
temp5 = data.frame()

# create progress bar
  
MS6 = rep(0, nrow(anovaTable6))
MS5 = rep(0, nrow(anovaTable5))
#MSCRD = rep(0, nrow(anovaTableCRD))
save.y6 = rep(0, length(y6))
save.y5 = rep(0, length(y6))

   
cat("gamma.ani = ", gamma.ani, "\n")
pb <- txtProgressBar(min = 0, max = nSim, style = 3)
counter = 0
# simulation

while (counter < nSim) {
    setTxtProgressBar(pb, counter + 1)
    
    run.eff = rnorm(nlevels(design.df$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
    cag.eff = rnorm(nlevels(design.df$Cag), mean = 0, sd = sqrt(gamma.cag * VC.resid))
    ani.eff = rnorm(nlevels(design.df$Ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
    trt.eff = c(1,2,3,4,5,6) #runif(nlevels(design.df$Trt), 0, 10)
    tag.eff = runif(nlevels(design.df$Tag), 0, 1)
    res.eff = rnorm(nrow(design.df), mean = 0, sd = sqrt(VC.resid))
    gm = 10
    
    
    y6 = gm + with(design6, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff
    
    save.y6 =  cbind(save.y6, y6)
    
    y5 = gm + with(design5, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff
    
    save.y5 =  cbind(save.y5, y5)
                
    yCRD = gm + with(designCRD, run.eff[Run] + ani.eff[Ani] + 
        tag.eff[Tag] + trt.eff[Trt]) + res.eff
                        


    aov.table6 = summaryAovTwoPhase(design.df = design6, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
      response = y6)
    
    aov.table5 = summaryAovTwoPhase(design.df = design5, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
      response = y5)
    
        
    #aov.table14
    
    #summary(aov(y14~ Tag + Trt+ Error(Run + Cag/Ani), data = design14))

    #aov.table16

    #summary(aov(y16~ Tag + Trt+ Error(Run + Cag/Ani), data = design16))
       
    
    #summary.lm(aov(y14~ Tag + Trt+ Error(Run + Cag/Ani), data = design14)$"Cag:Ani")
    #summary.lm(aov(y16~ Tag + Trt+ Error(Run + Cag/Ani), data = design16)$"Cag:Ani")

    #MSCRD = cbind(MSCRD, as.numeric(aov.tableCRD$A[, "MS"]))

    
    MS6 = cbind(MS6, as.numeric(aov.table6$A[, "MS"]))
    MS5 = cbind(MS5, as.numeric(aov.table5$A[, "MS"]))
    
    counter = counter + 1
}
close(pb)
proc.time() - old
 
MS6 = MS6[, -1]
MS5 = MS5[, -1]
  
 save.y6 =   save.y6[, -1]
  save.y5 =   save.y5[, -1]
    
            
mean(MS6[11, ]/(6 * 2600/3169))
mean(MS5[13, ]/(6 * 13160/15859))
       
mean(MS6[14,])
mean(MS5[16,])

mean(MS6[11,])
mean(MS5[13,])


mean((MS6[11, ] - MS6[14,])/2)
mean((MS5[13, ] - MS5[16,])/2)

sd((MS6[11, ] - MS6[14,])/2)
sd((MS5[13, ] - MS5[16,])/2)

mean(MS6[10, ]/MS6[11, ])
mean(MS5[12, ]/MS5[13, ])


sum(1 - pf(MS6[10, ]/MS6[11, ], 5, 6) < 0.05)/ nSim

sum(1 - pf(MS5[12, ]/MS5[13, ], 5, 5) < 0.05)/ nSim

sum(pf(MS6[10, ]/MS6[11, ], 5, 6)< 0.05/ nSim)

sum(p.adjust(1 - pf(MS5[12, ]/MS5[13, ], 5, 5), "fdr") < 0.05)/nSim

        
trtVar = sapply(MS6[11, ], function(x) sqrt(2 * x/(6 * contr6.e)))
trt.Est6 = apply(save.y6, 2, function(x) 2 * as.numeric(trtEff6 %*% x) %*% pairwiseContr(design6))
p.value6 = pt(-abs(trt.Est6/trtVar), df = 6) * 2
power6 = apply(p.value6, 1, function(x) sum(x < 0.05)/length(x))

effPowerMat6 = matrix(0, 6, 6)

for( i in 1:length(contr6.e)){   
  effPowerMat6[states[i,1], states[i,2]] = power6[i]
  effPowerMat6[states[i,2], states[i,1]] = contr6.e[i]
}


          [,1]      [,2]      [,3]      [,4]      [,5]   [,6]
[1,] 0 & 0.0932 & 0.2301 & 0.4466 & 0.6446 & 0.8664
[2,] 0.7738 & 0 & 0.0894 & 0.1917 & 0.4470 & 0.6040
[3,] 0.8759 & 0.8104 & 0 & 0.0989 & 0.2180 & 0.4538
[4,] 0.8744 & 0.7276 & 0.9075 & 0 & 0.0840 & 0.2534
[5,] 0.8104 & 0.8759 & 0.7738 & 0.7504 & 0 & 0.0860
[6,] 0.9075 & 0.7504 & 0.8744 & 0.9600 & 0.7276 & 0


trtVar = sapply(MS5[13, ], function(x) sqrt(2 * x/(6 * contr5.e)))
trt.Est5 = apply(save.y5, 2, function(x) 2 * as.numeric(trtEff5 %*% x) %*% pairwiseContr(design5))
p.value5 = pt(-abs(trt.Est5/trtVar), df = 5) * 2
power5 = apply(p.value5, 1, function(x) sum(x < 0.05)/length(x))

effPowerMat5 = matrix(0, 6, 6)

for( i in 1:length(contr5.e)){   
  effPowerMat5[states[i,1], states[i,2]] = power5[i]
  effPowerMat5[states[i,2], states[i,1]] = contr5.e[i]
}


[1,] 0 & 0.0871 & 0.2128 & 0.4175 & 0.5970 & 0.8363 \\
[2,] 0.8257 & 0 & 0.0809 & 0.1944 & 0.3900 & 0.6167 \\
[3,] 0.8751 & 0.7805 & 0 & 0.0962 & 0.2146 & 0.4171 \\
[4,] 0.8711 & 0.7701 & 0.9127 & 0 & 0.0930 & 0.2184 \\
[5,] 0.7805 & 0.7871 & 0.8257 & 0.8175 & 0 & 0.0864 \\
[6,] 0.9127 & 0.8175 & 0.8711 & 0.8710 & 0.7701 & 0 \\







 
#########################################################################################

test = 
function(){ 
 run.eff = rnorm(nlevels(design.df$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
    cag.eff = rnorm(nlevels(design.df$Cag), mean = 0, sd = sqrt(gamma.cag * VC.resid))
    ani.eff = rnorm(nlevels(design.df$Ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
    trt.eff = c(1,1,1,6,6,6) #runif(nlevels(design.df$Trt), 0, 10)
    tag.eff = runif(nlevels(design.df$Tag), 0, 1)
    res.eff = rnorm(nrow(design.df), mean = 0, sd = sqrt(VC.resid))
    gm = 10
    
    
    y6 = gm + with(design6, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff
       
    y5 = gm + with(design5, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff
                       
    aov.table6 = summaryAovTwoPhase(design.df = design6, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
      response = y6)
    
    aov.table5 = summaryAovTwoPhase(design.df = design5, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
      response = y5)
    
    return(c(as.numeric(aov.table6$A[, "MS"]),as.numeric(aov.table5$A[, "MS"]))) 
}

old = proc.time()
temp = foreach(i=1:nSim, .combine=rbind) %dopar% test()
proc.time() - old 

#######################################################################################
    y6 = gm + with(design6, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff
    
    y5 = gm + with(design5, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff
    
    aov.table6 = summaryAovTwoPhase(design.df = design6, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
      response = y6)
    
    aov.table5 = summaryAovTwoPhase(design.df = design5, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
      response = y5)

gamma.run =  10
gamma.cag =  7
gamma.ani =  5 #  c(10^((-4:4)/2))
 
 real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag *  VC.resid), (gamma.run * VC.resid))
  
 MS6.ChiSqu = t(suppressWarnings(apply(aov.table6$A, 1, 
          function(x) sum(fracToNum(x[2:5]) * real.VC) * rchisq(100000, fracToNum(x[1]))/fracToNum(x[1]))))

 MS5.ChiSqu = t(suppressWarnings(apply(aov.table5$A, 1, 
          function(x) sum(fracToNum(x[2:5]) * real.VC) * rchisq(100000, fracToNum(x[1]))/fracToNum(x[1]))))
 
 varDist( MS6ChiSqu,  MS5.ChiSqu, sample( 1:10000, 10000))

 
mean(MSCRD[9, ])         
mean(MS6[11, ])
mean(MS5[13, ])

mean(MSCRD[12, ])         
mean(MS6[14,])
mean(MS5[16,])

mean((MSCRD[9, ] -   MSCRD[12,])/2)
mean((MS6[11, ] -   MS6[14,])/2)
mean((MS5[13, ] -   MS5[16,])/2)


############################################################################################
 real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.cag *  VC.resid), (gamma.run * VC.resid))      


VCDist = function(MS6, MS5, interval){

  plot1 = qplot(MS6[11, interval], geom="histogram", binwidth=.5, xlab = "VC of measurement error with 6 DF")
  plot2 = qplot(MS5[13, interval], geom="histogram", binwidth=.5, xlab = "VC of measurement error with 5 DF")
  
  plot4 = qplot((MS6[11,interval]- MS6[14,interval])/2, geom="histogram", binwidth=.5, xlab = "VC of animals within cages with 6 DF")
  plot5 = qplot((MS5[13,interval] -  MS5[16,interval])/2, geom="histogram", binwidth=.5, xlab = "VC of animals within cages with 5 DF")
  
  plot7 = qplot((MS6[7,interval]- MS6[11,interval])/12, 
                    geom="histogram", binwidth=.5, xlab = "VC of cages with 6 DF")
  plot8 = qplot((MS5[9,interval ]-MS5[13, interval])/11, 
                geom="histogram", binwidth=.5, xlab = "VC of cages with 5 DF")
  
  plot9 = qplot((MS6[5,interval]- MS6[14,interval])/4, 
                    geom="histogram", binwidth=.5, xlab = "VC of runs with 6 DF")
  plot10 = qplot((MS5[6,interval ]-MS5[13, interval])/4, 
                geom="histogram", binwidth=.5, xlab = "VC of runs with 5 DF")

  grid.arrange(plot1, plot2,  plot4, plot5,  plot7, plot8,  plot9, plot10, ncol=2)
 

  result = data.frame("TRUE.value" = c(VC.resid, (gamma.ani * VC.resid),  (gamma.cag *  VC.resid), (gamma.run * VC.resid)),
  DF6mean  = c(mean(MS6[14, interval]),  mean((MS6[11,interval]- MS6[14,interval])/2), 
          mean((MS6[7,interval]- MS6[11,interval])/12),
          mean((MS6[5,interval]- MS6[14,interval])/4)),
  
  DF5mean  =c(mean(MS5[16, interval]), mean((MS5[13,interval] -  MS5[16,interval])/2), 
          mean((MS5[9,interval ]-MS5[13, interval])/11),
           mean((MS5[6,interval ]-MS5[16, interval])/4)))
 
  
  rownames(result) = c("Measurement error", "Animals within cages", "Cages", "Runs")
  
  return(result)
}

VCDist(MS6, MS5, sample(1:nSim, 1000))

getwd()
pdf(file = "VC6trt.pdf")
VCDist(MS6, MS5, sample(1:nSim, 100000))
dev.off()

################################################################################

dat <- data.frame(xx = c((MS6[11,]- MS6[14,])/2,(MS5[13,] -  MS5[16,])/2), Design = rep(c("6DF", "5DF"),each = nSim))

pdf(file = "VC6trtDist.pdf")
qplot(xx, data=dat, geom="density", fill=Design, alpha=I(.5), 
   main="Distribution of animal VC estimated from two designs", 
   xlab= expression(sigma[A]^2), 
   ylab="Density")
dev.off()

 mean((MS6[11,]- MS6[14,])/2)
sd((MS6[11,]- MS6[14,])/2)
 mean((MS5[13,] -  MS5[16,])/2)
sd((MS5[13,] -  MS5[16,])/2)


#################################################################################

trtDist = function(MS6, MS5, interval){

  plot1 = qplot(MS6[11,interval]/(6 * 2600/3169), geom="histogram", binwidth=.1, xlab = "Treatment variance with 6 DF")
  plot2 = qplot(MS5[13,interval]/(6 * 13160/15859), geom="histogram", binwidth=.1, xlab = "Treatment variance with 5 DF")
  
   
  plot7 = qplot(MS6[10, interval]/MS6[11,interval ], 
                    geom="histogram", binwidth=.5, xlab = "F-ratio with 6 DF")
  plot8 = qplot(MS5[12,interval ]/MS5[13, interval], 
                geom="histogram", binwidth=.5, xlab = "F-ratio with 5 DF")
  
   plot4 = qplot(p.adjust(1 - pf(MS6[10, interval]/MS6[11, interval], 5, 6), "fdr"), geom="histogram", binwidth=.05, xlab = "Adjusted p-value with 6 DF")
   
  plot5 = qplot(p.adjust(1 - pf(MS5[12, interval]/MS5[13, interval], 5, 5), "fdr"), geom="histogram", binwidth=.05, xlab = "Adjusted p-value with 5 DF")   

  grid.arrange(plot1, plot2, plot7, plot8, plot4, plot5,  ncol=2)
 

  result = data.frame(
  DF6mean  = c(mean(MS6[11,interval]/(6 * 2600/3169)), 
          mean(MS6[10, interval]/MS6[11,interval ]),
          sum(p.adjust(1 - pf(MS6[10, interval]/MS6[11, interval], 5, 6), "fdr") < 0.05)/ length(interval)),
  
  DF5mean  =c(mean(MS5[13,interval]/(6 * 13160/15859)), 
          mean(MS5[12,interval ]/MS5[13, interval]),
          sum(p.adjust(1 - pf(MS5[12, interval]/MS5[13, interval], 5, 5), "fdr") < 0.05)/length(interval)))
 
  
  rownames(result) = c("Treatment variance", "F-ratio", "Power of test")
  
  return(result)
}


trtDist(MS6, MS5, interval = sample(1:nSim, 10000))

getwd()
pdf(file = "trt6trt.pdf")
trtDist(MS6, MS5, sample(1:nSim, 100000))
dev.off()


#Pairwise comparison

contr.e = contrCompare(design6, pairwiseContr(design6))

contrastMat6 = matrix(0, ncol= 6, nrow = 6)

for( i in 1:length(contr.e)){   
  contrastMat6[states[i,2], states[i,1]] = contr.e[i]
}

trtVar = rowMeans(sapply(MS6[11,], function(x) sqrt(2*x/(6* contr.e)))) 

for( i in 1:length(trtVar)){   
  contrastMat6[states[i,1], states[i,2]] = trtVar[i]
}

contrastMat6

 
 trt.Est6 = apply(save.y6, 2, function(x) 2*tapply(x, design6$Trt, mean) %*% pairwiseContr(design6))

 p.value6 =  apply(trt.Est6, 2, function(x)  pt(-abs(x/trtVar), df = 6) * 2)

 power6 = apply(p.value6, 1, function(x) sum(p.adjust(x, "fdr")<0.05)/length(x)) 

 #power6 = apply(p.value6, 1, function(x) sum(x<0.05)/length(x)) 

 powerMat6 = matrix(0, ncol= 6, nrow = 6)

for( i in 1:length(contr.e)){   
  powerMat6[states[i,2], states[i,1]] = power6[i]
}
  
  trt.Est6 = apply(save.y6, 2, function(x) 2*tapply(x, design6$Trt, mean) %*% pairwiseContr(design6) * contr.e)

  trtEst6 = rowMeans(trt.Est6)
  
for( i in 1:length(trtVar)){   
  powerMat6[states[i,1], states[i,2]] = trtEst6 [i]
}

 powerMat6
 
 
contr.e = contrCompare(design5, pairwiseContr(design5))

contrastMat5 = matrix(0, ncol= 6, nrow = 6)

for( i in 1:length(contr.e)){   
contrastMat5[states[i,2], states[i,1]] = contr.e[i]
}

trtVar = rowMeans(sapply(MS5[13,], function(x) sqrt(2*x/(6* contr.e)))) 


for( i in 1:length(trtVar)){   
  contrastMat5[states[i,1], states[i,2]] = trtVar[i]
}

contrastMat5


 trt.Est5 = apply(save.y5, 2, function(x) 2*tapply(x, design5$Trt, mean) %*% pairwiseContr(design5))

 p.value5 =  apply(trt.Est5, 2, function(x)  pt(-abs(x/trtVar), df = 5) * 2)

 power5 = apply(p.value5, 1, function(x) sum(p.adjust(x, "fdr")<0.05)/length(x)) 
 #power5 = apply(p.value5, 1, function(x) sum(x<0.05)/length(x)) 

 powerMat5 = matrix(0, ncol= 6, nrow = 6)

for( i in 1:length(contr.e)){   
  powerMat5[states[i,2], states[i,1]] = power5[i]
}
  
  trt.Est5 = apply(save.y5, 2, function(x) 2*tapply(x, design5$Trt, mean) %*% pairwiseContr(design5) * contr.e)

  trtEst5 = rowMeans(trt.Est5)
  
for( i in 1:length(trtVar)){   
  powerMat5[states[i,1], states[i,2]] = trtEst5[i]
}

 powerMat5
 powerMat6

#Pick the best 3 pairwise contrasts 

contr1 = c(1, 0, 0,0,0,-1)
contr2 = c(0, 0,1, -1,0,0)
contr3 = c( 0,1,0,0,-1, 0)
contr4 = c(1, 0,-1,-1,0,1)
contr5 = c(0.5,-1, 0.5,0.5,-1,0.5)
  
contr = cbind(contr1, contr2, contr3, contr4, contr5)/2  
contr

contrCompare(design5, contr)



sqrt(2*as.numeric(MS6[11,])/(6* contrCompare(design6, contr)))


##############################################################################################
      
contr1 = c(1, 1, 1,-1,-1,-1)
contr2 =c(0.5, 0.5, -1,0.5,0.5,-1)
contr3 = c( 1, -1,0, 1, -1, 0)
contr4 = contr1 * contr2
contr5 = contr1 * contr3

contr = cbind(contr1, contr2, contr3, contr4, contr5)/2  
contr


     contr61 = contr1[design6$Trt]/2
    contr62 = contr2[design6$Trt]/2
    contr63 = contr3[design6$Trt]/2
    contr64 = contr4[design6$Trt]/2
    contr65 = contr5[design6$Trt]/2


anovaTable6 = summaryAovTwoPhase(design.df = design6, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
    trt.contr =  list(Tag = NA,
    Trt = list(contr1 = contr61, contr2 = contr62, 
    contr3 = contr63, contr4 = contr64, contr5 = contr65)),
    contr.matrix = TRUE,  response = y6)$A

    contr51 = contr1[design5$Trt]/2
    contr52 = contr2[design5$Trt]/2
    contr53 = contr3[design5$Trt]/2
    contr54 = contr4[design5$Trt]/2
    contr55 = contr5[design5$Trt]/2

anovaTable5 = summaryAovTwoPhase(design.df = design5, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", trt.contr =  list(Tag = NA, 
        Trt = list(contr1 = contr51, contr2 = contr52, contr3 = contr53, 
        contr4 = contr54, contr5 = contr55)),
    contr.matrix = TRUE, 
      response = y5)$A

VC.resid = 1

gamma.run =  10
gamma.cag =  5
gamma.ani =  5 #  c(10^((-4:4)/2))

nSim = 10000


# create progress bar
  
MS6 = rep(0, nrow(anovaTable6))
MS5 = rep(0, nrow(anovaTable5))
   
cat("gamma.ani = ", gamma.ani, "\n")
pb <- txtProgressBar(min = 0, max = nSim, style = 3)
counter = 0
# simulation

while (counter < nSim) {
    setTxtProgressBar(pb, counter + 1)
    
    run.eff = rnorm(nlevels(design.df$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
    cag.eff = rnorm(nlevels(design.df$Cag), mean = 0, sd = sqrt(gamma.cag * VC.resid))
    ani.eff = rnorm(nlevels(design.df$Ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
    trt.eff = c(1, 6, 3.5, 1, 6, 3.5) #runif(nlevels(design.df$Trt), 0, 10)
    tag.eff = runif(nlevels(design.df$Tag), 0, 1)
    res.eff = rnorm(nrow(design.df), mean = 0, sd = sqrt(VC.resid))
    gm = 10
    
    
    y6 = gm + with(design6, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff
    
    y5 = gm + with(design5, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff

    contr61 = contr1[design6$Trt]/2
    contr62 = contr2[design6$Trt]/2
    contr63 = contr3[design6$Trt]/2
    contr64 = contr4[design6$Trt]/2
    contr65 = contr5[design6$Trt]/2
    
    aov.table6 = summaryAovTwoPhase(design.df = design6, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
    trt.contr =  list(Tag = NA, Trt = list(contr1 = contr61, contr2 = contr62, 
                                            contr3 = contr63, contr4 = contr64, contr5 = contr65)),
    contr.matrix = TRUE,  response = y6)
        
    
    contr51 = contr1[design5$Trt]/2
    contr52 = contr2[design5$Trt]/2
    contr53 = contr3[design5$Trt]/2
    contr54 = contr4[design5$Trt]/2
    contr55 = contr5[design5$Trt]/2

    aov.table5 = summaryAovTwoPhase(design.df = design5, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", trt.contr =  list(Tag = NA, 
        Trt = list(contr1 = contr51, contr2 = contr52, contr3 = contr53, 
        contr4 = contr54, contr5 = contr55)),
    contr.matrix = TRUE, 
      response = y5)
      
        #aov.table14
    
    #summary(aov(y14~ Tag + Trt+ Error(Run + Cag/Ani), data = design14))

    #aov.table16

    #summary(aov(y16~ Tag + Trt+ Error(Run + Cag/Ani), data = design16))
       
    
    #summary.lm(aov(y14~ Tag + Trt+ Error(Run + Cag/Ani), data = design14)$"Cag:Ani")
    #summary.lm(aov(y16~ Tag + Trt+ Error(Run + Cag/Ani), data = design16)$"Cag:Ani")

   # MSCRD = cbind(MSCRD, as.numeric(aov.tableCRD$A[, "MS"]))

    MS6 = cbind(MS6, as.numeric(aov.table6$A[, "MS"]))
    MS5 = cbind(MS5, as.numeric(aov.table5$A[, "MS"]))
    
    counter = counter + 1
}
close(pb)

 
MS6 = MS6[, -1]
MS5 = MS5[, -1]
  
          
mean(MS6[11, ])
mean(MS5[13, ])
       
mean(MS6[14,])
mean(MS5[16,])

mean((MS6[11, ] - MS6[14,])/2)
mean((MS5[13, ] - MS5[16,])/2)



sum(1 - pf(MS6[10, ]/MS6[11, ], 5, 6) < 0.05)/ nSim

sum(1 - pf(MS5[12, ]/MS5[13, ], 5, 5) < 0.05)/ nSim

sum(p.adjust(1 - pf(MS6[10, ]/MS6[11, ], 5, 6), "fdr") < 0.05)/ nSim

sum(p.adjust(1 - pf(MS5[12, ]/MS5[13, ], 5, 5), "fdr") < 0.05)/nSim


sum(p.adjust(1 - pf(MS6[14,]/MS6[17,], 1, 6), "fdr") < 0.05)/nSim


sum(p.adjust(1 -pf( MS5[17,]/MS5[20, ], 1, 5), "fdr") < 0.05)/nSim

MS5[17,]/MS5[20, ]

############################################################################################
#Compute the SED

   run.eff = rnorm(nlevels(design.df$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
    cag.eff = rnorm(nlevels(design.df$Cag), mean = 0, sd = sqrt(gamma.cag * VC.resid))
    ani.eff = rnorm(nlevels(design.df$Ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
    trt.eff = c(1,1,1,6,6,6) #runif(nlevels(design.df$Trt), 0, 10)
    tag.eff = runif(nlevels(design.df$Tag), 0, 1)
    res.eff = rnorm(nrow(design.df), mean = 0, sd = sqrt(VC.resid))
    gm = 10
    
    
    y6 = gm + with(design6, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff
    
    y5 = gm + with(design5, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff

    (aov.table6 = summaryAovTwoPhase(design.df = design6, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
      response = y6)$A)
   
 aov6.fit = aov(y6~  Tag + Trt + Error(Run + Cag/Ani), design6)
  summary(aov6.fit)

  summary.lm(aov6.fit$"Cag:Ani")

  
 pt(-abs(summary.lm(aov6.fit$"Cag:Ani")$coef[,"t value"][-1]), df = 6) * 2
 
 
  
contr1 = c(-1, 1, 0,0,0,0)
contr2 =c(-1, 0, 1,0,0,0)
contr3 = c(- 1, 0,0, 1, 0, 0)
contr4 =c( -1, 0,0, 0, 1, 0)
contr5 = c(- 1, 0,0, 0, 0, 1)


contr = cbind(contr1, contr2, contr3, contr4, contr5)/2  
contr

contrCompare(design6, contr)


sqrt(2*as.numeric(aov.table6[11,"MS"])/(6* contrCompare1(design6, contr))  )
  
  
 as.numeric(trtEff6 %*% y6  ) %*% contr   *2
  
  
apply(contr, 2, function(x) sum(-x * tapply(y6, design6$Trt, mean)))  
  
  contrCompare(design6, basicContr(design6)[,-6])
  
  
   (aov.table5 = summaryAovTwoPhase(design.df = design6, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
      response = y5)$A)
   
 aov5.fit = aov(y5~  Tag + Trt + Error(Run + Cag/Ani), design5)
  summary(aov5.fit)

  summary.lm(aov5.fit$"Cag:Ani")
 
  
##################################################################################################  
#Testing the Power
    groups = length(trt.eff) 
    between.var = var(trt.eff)
    within.var = 1 + 2*2
    
    sig.level = 0.05
    lambda <- (groups  -1 ) *  6 * 2600/3169 * (between.var/within.var)
      lower.df = 6 
 
     #lambda <- (groups - 1) *  MS6[10,] /(MS6[11, ]/(2600/3169))
     
    pf(qf(sig.level, groups - 1, lower.df , lower.tail = FALSE), 
            groups - 1, lower.df , lambda, lower.tail = FALSE)
            
           
      lower.df = 5 
      lambda <- (groups - 1)  * 6  *  13160/15859 * (between.var/within.var)
      
     pf(qf(sig.level, groups - 1, lower.df , lower.tail = FALSE), 
            groups - 1, lower.df , lambda, lower.tail = FALSE)
    
   
     sum(p.adjust(pf(MS5[12,]/MS5[13,], 5, 5, lower.tail = FALSE), "fdr") < 0.05)/nSim
      
    (n - 1) * groups
    
  plot1 = qplot(MS6[11,interval]/(6 * 2600/3169), geom="histogram", binwidth=.1, xlab = "Treatment variance with 6 DF")
  plot2 = qplot(MS5[13,interval]/(6 * 13160/15859), geom="histogram", binwidth=.1, xlab = "Treatment variance with 5 DF")
 
 
 
 pt(qt(0.05/2, 6, lower.tail = FALSE), 6, ncp = sqrt(n/tsample) * 
            delta/sd, lower.tail = FALSE)
 
##############################################################################################
      
contr1 = c(1, 1, 1,-1,-1,-1)
contr2 =c(0.5, 0.5, -1,0.5,0.5,-1)
contr3 = c( 1, -1,0, 1, -1, 0)
contr4 = contr1 * contr2
contr5 = contr1 * contr3

contr = cbind(contr1, contr2, contr3, contr4, contr5)/2  
contr


     contr61 = contr1[design6$Trt]/2
    contr62 = contr2[design6$Trt]/2
    contr63 = contr3[design6$Trt]/2
    contr64 = contr4[design6$Trt]/2
    contr65 = contr5[design6$Trt]/2


anovaTable6 = summaryAovTwoPhase(design.df = design6, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
    trt.contr =  list(Tag = NA,
    Trt = list(contr1 = contr61, contr2 = contr62, 
    contr3 = contr63, contr4 = contr64, contr5 = contr65)),
    contr.matrix = FALSE,  response = y6)$A

anovaTable6 = summaryAovTwoPhase(design.df = design6, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
    trt.contr =  list(Tag = NA,
    Trt = list(contr1 = contr61)),
    contr.matrix = FALSE,  response = y6)


    contr51 = contr1[design5$Trt]/2
    contr52 = contr2[design5$Trt]/2
    contr53 = contr3[design5$Trt]/2
    contr54 = contr4[design5$Trt]/2
    contr55 = contr5[design5$Trt]/2

anovaTable5 = summaryAovTwoPhase(design.df = design5, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", trt.contr =  list(Tag = NA, 
        Trt = list(contr1 = contr51, contr2 = contr52, contr3 = contr53, 
        contr4 = contr54, contr5 = contr55)),
    contr.matrix = FALSE, 
      response = y5)$A

VC.resid = 1

gamma.run =  10
gamma.cag =  5
gamma.ani =  2 #  c(10^((-4:4)/2))

nSim = 10000


# create progress bar
  
MS6 = rep(0, nrow(anovaTable6))
MS5 = rep(0, nrow(anovaTable5))
   
cat("gamma.ani = ", gamma.ani, "\n")
pb <- txtProgressBar(min = 0, max = nSim, style = 3)
counter = 0
# simulation

while (counter < nSim) {
    setTxtProgressBar(pb, counter + 1)
    
    run.eff = rnorm(nlevels(design.df$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
    cag.eff = rnorm(nlevels(design.df$Cag), mean = 0, sd = sqrt(gamma.cag * VC.resid))
    ani.eff = rnorm(nlevels(design.df$Ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
    trt.eff = c(1, 6, 3.5, 1, 6, 3.5) #runif(nlevels(design.df$Trt), 0, 10)
    tag.eff = runif(nlevels(design.df$Tag), 0, 1)
    res.eff = rnorm(nrow(design.df), mean = 0, sd = sqrt(VC.resid))
    gm = 10
    
    
    y6 = gm + with(design6, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff
    
    y5 = gm + with(design5, run.eff[Run] + ani.eff[ani] + cag.eff[Cag] + 
            tag.eff[Tag] + trt.eff[Trt]) + res.eff

    contr61 = contr1[design6$Trt]/2
    contr62 = contr2[design6$Trt]/2
    contr63 = contr3[design6$Trt]/2
    contr64 = contr4[design6$Trt]/2
    contr65 = contr5[design6$Trt]/2
    
    aov.table6 = summaryAovTwoPhase(design.df = design6, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", 
    trt.contr =  list(Tag = NA, Trt = list(contr1 = contr61, contr2 = contr62, 
                                            contr3 = contr63, contr4 = contr64, contr5 = contr65)),
    contr.matrix = FALSE,  response = y6)
        
    
    contr51 = contr1[design5$Trt]/2
    contr52 = contr2[design5$Trt]/2
    contr53 = contr3[design5$Trt]/2
    contr54 = contr4[design5$Trt]/2
    contr55 = contr5[design5$Trt]/2

    aov.table5 = summaryAovTwoPhase(design.df = design5, blk.str2 = "Run", 
    blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", trt.contr =  list(Tag = NA, 
        Trt = list(contr1 = contr51, contr2 = contr52, contr3 = contr53, 
        contr4 = contr54, contr5 = contr55)),
    contr.matrix = FALSE, 
      response = y5)
      
        #aov.table14
    
    #summary(aov(y14~ Tag + Trt+ Error(Run + Cag/Ani), data = design14))

    #aov.table16

    #summary(aov(y16~ Tag + Trt+ Error(Run + Cag/Ani), data = design16))
       
    
    #summary.lm(aov(y14~ Tag + Trt+ Error(Run + Cag/Ani), data = design14)$"Cag:Ani")
    #summary.lm(aov(y16~ Tag + Trt+ Error(Run + Cag/Ani), data = design16)$"Cag:Ani")

   # MSCRD = cbind(MSCRD, as.numeric(aov.tableCRD$A[, "MS"]))

    MS6 = cbind(MS6, as.numeric(aov.table6$A[, "MS"]))
    MS5 = cbind(MS5, as.numeric(aov.table5$A[, "MS"]))
    
    counter = counter + 1
}
close(pb)

 
MS6 = MS6[, -1]
MS5 = MS5[, -1]
  
          
mean(MS6[11, ])
mean(MS5[13, ])
       
mean(MS6[14,])
mean(MS5[16,])

mean((MS6[11, ] - MS6[14,])/2)
mean((MS5[13, ] - MS5[16,])/2)



sum(1 - pf(MS6[10, ]/MS6[11, ], 5, 6) < 0.05)/ nSim

sum(1 - pf(MS5[12, ]/MS5[13, ], 5, 5) < 0.05)/ nSim

sum(p.adjust(1 - pf(MS6[10, ]/MS6[11, ], 5, 6), "fdr") < 0.05)/ nSim

sum(p.adjust(1 - pf(MS5[12, ]/MS5[13, ], 5, 5), "fdr") < 0.05)/nSim



sum(p.adjust(1 - pf(MS6[14,]/MS6[17,], 1, 6), "fdr") < 0.05)/nSim


sum(p.adjust(1 -pf( MS5[17,]/MS5[20, ], 1, 5), "fdr") < 0.05)/nSim

MS5[17,]/MS5[20, ]















 design6 = cbind( design6, Sam = 1:36)
 summaryAovOnePhase(design.df = design6, blk.str = "Sam", trt.str = "Tag + Trt", trt.contr =  list(Tag = NA,
     Trt = list(contr1 = contr61, contr2 = contr62,
     contr3 = contr63, contr4 = contr64, contr5 = contr65)),     contr.matrix = FALSE,  response = y6) 



