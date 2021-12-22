n = 16
nPlot = 4

nCag = 2
nAni = 4

nTrt = 2

phase1DesignEX1 <- local({
  Cag = as.factor(rep(1:nCag, each = 2))
  Ani = as.factor(LETTERS[1:nAni])
  Trt = as.factor(letters[1:nTrt])
  data.frame(Cag, Ani,Trt)
})
phase1DesignEX1

summary.aov.onePhase(phase1DesignEX1,  blk.str = "Cag/Ani",  trt.str = "Trt")

#Parameter's of block
nBlk = n/nPlot

Zb = matrix(0, ncol = nBlk, nrow = n)
Z.des = rep(1:nBlk,  each = nPlot)
Zb[cbind(1:n, Z.des)] <- 1

Zt = matrix(0, ncol = nPlot, nrow = n)
tag.des = rep(1:nPlot,  time = nBlk)
Zt[cbind(1:n, tag.des)] <- 1

Pb = projMat(Zb)
Pb1 = projMat(Zt)

betRun = Pb - mK(n)
betTag = Pb1 - mK(n)
withBlock = (mI(n) - Pb) %*% (mI(n) - Pb1)

blk.proj = withBlock


#Parameter's of block structure of Phase 1
nZ1 = nAni
Z1.rep =  n/nZ1

Z1.des = phase1DesignEX1$Ani
Z1.mat = matrix(0, ncol = nZ1, nrow= nZ1*Z1.rep)
Z1.mat[cbind(1:( nZ1*Z1.rep), Z1.des)] = 1

#Parameter's of treatment
trt.rep = n/nTrt

trt.des = as.numeric(phase1DesignEX1$Trt[match(LETTERS[Z1.des], phase1DesignEX1$Ani)])
X.trt = matrix(0, ncol = nTrt, nrow= nTrt*trt.rep )
X.trt[cbind(1:( nTrt*trt.rep ), trt.des)] = 1

##################################################################################
#Function that swap any two random row indexes
swap<-function(ind){

  idx <- seq(2, length(ind)-2)
  changepoints <- sample(idx, size=2, replace=FALSE)
  tmp <- ind[changepoints[1]]
  ind[changepoints[1]] <- ind[changepoints[2]]
  ind[changepoints[2]] <- tmp
  return(ind)
 }

#Function to calculate the average efficiency factors
obj.fun.phase2=function(ind){
  Rep = Z1.rep
  newZ1 = Z1.mat[ind[-length(ind)],]
  info.mat =  t(newZ1) %*% blk.proj %*% newZ1
  e.va = eigen(info.mat)$va
  can.eff = e.va[-which(e.va<1e-7)]/Rep

  aveEff = 1/mean(1/can.eff) * 100
  return(aveEff)
}

################################################################################
# test run SANN this step need to repeat Many time to obtain a better results
# with different temperature

y2 = 0.1

init =  sample(1:n)
init = c(init, init[1])
res <- optim(init, obj.fun.phase2, swap, method="SANN",
    control =
      list(maxit = 10000, temp = y2 , tmax = 100, trace = TRUE,
                                                      REPORT = 1, fnscale = -1))

newind = res$par[-length(res$par)]

  #check the average efficiency factor of Phase 1 block assign to Phase 2 block]
test(X.trt = Z1.mat[newind,],  blk.proj, Rep=Z1.rep)$ave.eff

new.Z1.mat=  Z1.mat[newind,]

#Fit to the design
colnames(new.Z1.mat) = sort(levels(interaction(LETTERS[1:nZ1])))
new.Z1.des = apply(new.Z1.mat, 1,  function(x)
colnames(new.Z1.mat)[which(as.logical(x))])
 new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))

design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

design.df = cbind(design.df,
            Ani = t(new.Z1.des ),
            Cag = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                          as.character(phase1DesignEX1$Ani)),]$Cag,
            Trt = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                          as.character(phase1DesignEX1$Ani)),]$Trt)

design.df

summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani", trt.str = "Trt + Tag")
summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Trt + Tag")

summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")

c(nTrt, nAni, nCag, nPlot, n/nPlot)
fileName = paste("design_", nTrt, "trt_", nAni, "ani_", nCag, "Cag_", nPlot, "tag_", n/nPlot, "run.Rdata", sep = "")

save(design.df, file =fileName)

################################################################################
# run SANN this step need to repeat Many time to obtain a better results
# with different temperature

nIter = 100
y2 = 0.1

bestInd = sample(1:n)
bestAveEff = 0.0

pb <- txtProgressBar(min = 0, max = nIter, style = 3)

for(i in 1:nIter){
  setTxtProgressBar(pb, i)
  init =  sample(1:n)
  init = c(init, init[1])
  res <- optim(init, obj.fun.phase2, swap, method="SANN",
    control =
      list(maxit = 500000, temp = y2 , tmax = 1000, trace = FALSE,
                                                      REPORT = 1, fnscale = -1))

  #tempList = 1 / log((((1:100000000)-1) %/% 10000)*10000 + exp(1))

  newind = res$par[-length(res$par)]

  #check the average efficiency factor of Phase 1 block assign to Phase 2 block]
  oldAveEff = test(X.trt = Z1.mat[newind,],  blk.proj, Rep=Z1.rep)$ave.eff
  if(bestAveEff < oldAveEff){
    bestAveEff = oldAveEff
    bestInd = newind
  }
}
close(pb)
