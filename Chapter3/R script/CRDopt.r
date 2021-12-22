
testCRD = function(nTrt, nAni, nTag, nRun)

nTrt = 3
bRep = 2
tRep = 2
n = nTrt * bRep * tRep
nPlot = 4
nBlk = n/nPlot

nAni = nTrt * bRep

c(nTrt, nAni, nPlot, n/nPlot)

phase1DesignEX1 <- local({
    Ani = as.factor(LETTERS[1:nAni])
    Trt = as.factor(letters[1:nTrt])
    data.frame(Ani, Trt)
})
phase1DesignEX1

summary.aov.onePhase(phase1DesignEX1, blk.str = "Ani", trt.str = "Trt")

# Parameter's of block
nBlk = n/nPlot

Zb = matrix(0, ncol = nBlk, nrow = n)
Z.des = rep(1:nBlk, each = nPlot)
Zb[cbind(1:n, Z.des)] <- 1

Zt = matrix(0, ncol = nPlot, nrow = n)
tag.des = rep(1:nPlot, time = nBlk)
Zt[cbind(1:n, tag.des)] <- 1

Pb = projMat(Zb)
Pb1 = projMat(Zt)

betRun = Pb - mK(n)
betTag = Pb1 - mK(n)
withBlock = (mI(n) - Pb) %*% (mI(n) - Pb1)

blk.proj = withBlock


# Parameter's of block structure of Phase 1
nZ1 = nrow(phase1DesignEX1)
Z1.rep = n/nZ1

Z1.des = as.factor(rep(1:nZ1, Z1.rep))
Z1.mat = matrix(0, ncol = nZ1, nrow = nZ1 * Z1.rep)
Z1.mat[cbind(1:(nZ1 * Z1.rep), Z1.des)] = 1


#Z1.mat = xxxx$b

# Parameter's of treatment
trt.rep = n/nTrt

trt.des = as.numeric(phase1DesignEX1$Trt[match(LETTERS[Z1.des], phase1DesignEX1$Ani)])
X.trt = matrix(0, ncol = nTrt, nrow = nTrt * trt.rep)
X.trt[cbind(1:(nTrt * trt.rep), trt.des)] = 1 


swap <- function(ind) {
    
    idx <- seq(2, length(ind) - 2)
    changepoints <- sample(idx, size = 2, replace = FALSE)
    tmp <- ind[changepoints[1]]
    ind[changepoints[1]] <- ind[changepoints[2]]
    ind[changepoints[2]] <- tmp
    return(ind)
}

obj.fun.new = function(ind) {
    Rep = Z1.rep
    newZ1 = Z1.mat[ind[-length(ind)], ]
    info.mat = t(newZ1) %*% blk.proj %*% newZ1
    e.va = eigen(info.mat)$va
    can.eff = e.va[-which(e.va < 1e-07)]/Rep
    
    aveEff = 1/mean(1/can.eff)
   
    Rep = trt.rep
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*% blk.proj %*% newX
    e.va = eigen(info.mat)$va
    can.eff = e.va[-which(e.va < 1e-07)]/Rep

    aveEff = aveEff * (3/4) + (1/mean(1/can.eff)) *(1/4)

    return(aveEff)
}

##################################################################################
#Function that swap any two random row indexes
swap <- function(ind) {
    
    idx <- seq(2, length(ind) - 2)
    changepoints <- sample(idx, size = 2, replace = FALSE)
    tmp <- ind[changepoints[1]]
    ind[changepoints[1]] <- ind[changepoints[2]]
    ind[changepoints[2]] <- tmp
    return(ind)
}

length(can.eff)

# Function to calculate the average efficiency factors
obj.fun.ms = function(ind) {
    newZ1 = Z1.mat[ind[-length(ind)], ]
    info.mat = t(newZ1) %*% blk.proj %*% newZ1
    ms1 = tr(info.mat)
    ms2 = tr(info.mat^2)
    
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*% blk.proj %*% newX

    result =  ms1/ms2 + tr(info.mat)/tr(info.mat^2)
    return(result)
} 

obj.fun.ms = function(ind) {
    newZ1 = Z1.mat[ind[-length(ind)], ]
    info.mat = t(newZ1) %*% blk.proj %*% newZ1
    ms1 = tr(info.mat)
    ms2 = tr(info.mat^2)
   
    return(ms1/ms2)
} 

# Function to calculate the average efficiency factors
obj.fun.phase2 = function(ind) {
    Rep = Z1.rep
    newZ1 = Z1.mat[ind[-length(ind)], ]
    info.mat = t(newZ1) %*% blk.proj %*% newZ1
    e.va = eigen(info.mat)$va
    can.eff = e.va[-which(e.va < 1e-07)]/Rep
    
    aveEff = 1/mean(1/can.eff) * 100
    return(aveEff)
} 

##################################################################################

obj.fun.new = function(ind) {
    Rep = Z1.rep
    newZ1 = Z1.mat[ind[-length(ind)], ]
    info.mat = t(newZ1) %*% blk.proj %*% newZ1
    e.va = eigen(info.mat)$va
    can.eff = e.va[-which(e.va < 1e-07)]/Rep
    
    aveEff = 1/mean(1/can.eff) * 100
   
    Rep = trt.rep
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*% blk.proj %*% newX
    e.va = eigen(info.mat)$va
    can.eff = e.va[-which(e.va < 1e-07)]/Rep

    aveEff = aveEff * (3/4) + (1/mean(1/can.eff) * 100) *(1/4)

    return(aveEff)
} 

obj.fun.new1 = function(ind) {
    Rep = Z1.rep
    newZ1 = Z1.mat[ind[-length(ind)], ]
    info.mat = t(newZ1) %*% blk.proj %*% newZ1
    e.va = eigen(info.mat)$va
    can.eff = e.va[-which(e.va < 1e-07)]/Rep
    
    aveEff = 1/mean(1/can.eff) * 100 * length(can.eff)
   
    Rep = trt.rep
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*% blk.proj %*% newX
    e.va = eigen(info.mat)$va
    can.eff = e.va[-which(e.va < 1e-07)]/Rep

    aveEff = aveEff * (1/2) + (1/mean(1/can.eff) * 100* length(can.eff)) *(1/2)

    return(aveEff)
} 

obj.fun.phase1= function(ind){
  Rep = trt.rep
  newX = X.trt[ind[-length(ind)],]
  newZ1 = Z1.mat[ind[-length(ind)],]
  PP = blk.proj %*% newZ1 %*% (mI(nZ1) - mK(nZ1)) %*%
     invInfMat(blk.proj, newZ1, mI(nZ1) - mK(nZ1)) %*%
      (mI(nZ1) - mK(nZ1)) %*% t(newZ1) %*% blk.proj

  info.mat =  t(newX) %*% PP %*% newX
  e.va = Re(eigen(info.mat)$va)
  can.eff = e.va[-which(e.va<1e-7)]/Rep

  aveEff = 1/mean(1/can.eff)
  
  return(aveEff)
}


obj.fun.combine=function(ind){
  Rep = trt.rep
  newX = X.trt[ind[-length(ind)],]
  newZ1 = new.Z1.mat[ind[-length(ind)],]

  PP = blk.proj %*% newZ1 %*% (mI(nZ1) - mK(nZ1)) %*%
     invInfMat(blk.proj, newZ1, mI(nZ1) - mK(nZ1)) %*%
      (mI(nZ1) - mK(nZ1)) %*% t(newZ1) %*% blk.proj

  info.mat =  t(newX) %*% PP %*% newX
  e.va = Re(eigen(info.mat)$va)
  can.eff = e.va[-which(e.va<1e-7)]/Rep

  aveEff = 1/mean(1/can.eff)
  Rep = Z1.rep
  newZ1 = Z1.mat[ind[-length(ind)],]
  info.mat =  t(newZ1) %*% blk.proj %*% newZ1
  e.va = Re(eigen(info.mat)$va)
  can.eff = e.va[-which(e.va<1e-7)]/Rep

  aveEff = 0.5 * aveEff + 0.5* 1/mean(1/can.eff)

  return(aveEff)
}

###############################################################################
#Phase 2 optimisation
#initialised a better temperture for SA
states = lapply(1:nrow(X.trt), function(x) t(cbind(x, (x+1):nrow(X.trt))))
states = states[-length(states)]

#combine all the states into a matrix
states= t(matrix(c(states, recursive=TRUE), nrow = 2))

cpp.states = states - 1
C.trt.mat = diag(ncol(Z1.mat))

Z1.mat = Z1.mat[sample(nrow(Z1.mat)),]

#Chain 0
xx =  cpp.test(Z1.mat, blk.proj, cpp.states, C.trt.mat, Z1.rep, 10000, 10^9, 0)

temp = xx$storeAveEff[-which(xx$storeAveEff==0)]
temp = diff(temp)
temp = temp[which(temp<0)]
(y2 = mean(temp)/log(0.5))

y2  = y2 * 100




temp = replicate(10000, obj.fun.new(swap(init)), simplify = "vector")
temp = diff(temp)
temp = temp[which(temp<0)]
y2 = mean(temp)/log(runif(1, min=0, max=1))
                                  
################################################################################
# run SANN this step need to repeat Many time to obtain a better results
# with different temperature

y2 = 0.1

init =  sample(1:n)
init = c(init, init[1])
res <- optim(init, obj.fun.ms, swap, method="SANN",
    control =
      list(maxit = 100000, temp = y2 , tmax = 1000, trace = TRUE,
                                                      REPORT = 10, fnscale = -1))

newind = res$par[-length(res$par)]


y2 = 1

init =  sample(1:n)
init = c(init, init[1])
res <- optim(init, obj.fun.phase2, swap, method="SANN",
    control =
      list(maxit = 100000, temp = y2 , tmax = 1000, trace = TRUE,
                                                      REPORT = 10, fnscale = -1))

newind = res$par[-length(res$par)]



#Initialise the temperture for SA
init = sample(1:n)
init = c(init, init[1])

temp = replicate(100000, obj.fun.new(swap(init)), simplify = "vector")
temp = diff(temp)
temp = temp[which(temp<0)]

y2 = mean(temp/log(runif(length(temp), min=0, max=1)))

y2 = mean(temp)/log(runif(1, min=0, max=1))

y2 = mean(temp)/log(0.5)

#y2 = 8

init = sample(1:n)
init = c(init, init[1])

Rprof("example.out")
res <- optim(init, obj.fun.new, swap, method="SANN",
    control =
      list(maxit = 1000000, temp = y2 , tmax = 1000, trace = TRUE,
                                                     REPORT = 100, fnscale = -1))
Rprof(NULL)
summaryRprof("example.out") 

newind = res$par[-length(res$par)]
test(X.trt = Z1.mat[newind,], mI(n) - Pb, Rep=Z1.rep)$ave.eff
test(X.trt = Z1.mat[newind,], mI(n) - Pb, Rep=Z1.rep)$can.eff

y2 = y2 / (1 + 0.25)

res <- optim(res$par, obj.fun.new, swap, method="SANN",
    control =
      list(maxit = 100000, temp = y2 , tmax = 1000, trace = TRUE,
                                                     REPORT = 10, fnscale = -1))

newind = res$par[-length(res$par)]
test(X.trt = Z1.mat[newind,],  blk.proj, Rep=Z1.rep)$ave.eff
test(X.trt = Z1.mat[newind,], mI(n) - Pb, Rep=Z1.rep)$ave.eff

test(X.trt = Z1.mat[newind,], mI(n) - Pb, Rep=Z1.rep)$can.eff

test(X.trt = X.trt[newind,], mI(n) - Pb, Rep=trt.rep)$ave.eff
test(X.trt = X.trt[newind,], blk.proj, Rep=trt.rep)$can.eff

y2 = 1

init =  sample(1:n)
init = c(init, init[1])
res <- optim(init, obj.fun.new1, swap, method="SANN",
    control =
      list(maxit = 100000, temp = y2 , tmax = 1000, trace = TRUE,
                                                     REPORT = 10, fnscale = -1))

res <- optim(res$par, obj.fun.new1, swap, method="SANN",
    control =
      list(maxit = 100000, temp = y2 , tmax = 1000, trace = TRUE,
                                                     REPORT = 10, fnscale = -1))

newind = res$par[-length(res$par)]

#check the average efficiency factor of Phase 1 block assign to Phase 2 block]
test(X.trt = Z1.mat[newind,],  blk.proj, Rep=Z1.rep)$ave.eff
test(X.trt = Z1.mat[newind,],  blk.proj, Rep=Z1.rep)$can.eff

test(X.trt = Z1.mat[newind,], mI(n) - Pb, Rep=Z1.rep)$ave.eff
test(X.trt = Z1.mat[newind,],   mI(n) - Pb, Rep=Z1.rep)$can.eff

test(X.trt = Z1.mat[bestInd,],  blk.proj, Rep=Z1.rep)$ave.eff
test(X.trt = Z1.mat[newind,],  blk.proj, Rep=Z1.rep)$can

PP = (mI(n) - Pb) %*% Z1.mat[newind,] %*% (mI(nZ1) - mK(nZ1)) %*%
     invInfMat((mI(n) - Pb), Z1.mat[newind,], mI(nZ1) - mK(nZ1)) %*%
      (mI(nZ1) - mK(nZ1)) %*% t(Z1.mat[newind,]) %*% (mI(n) - Pb)

#check the average efficiency factor of traetments
test(X.trt = X.trt[newind,], blk.proj = PP, Rep = trt.rep )

(y2 = y2/(1 + 0.25)) #Cooling the temp for next chain

new.Z1.mat=  Z1.mat[newind,]

##################################################################################
nIter = 10
y2 = 0.1

indList = numeric(n)

pb <- txtProgressBar(min = 0, max = nIter, style = 3)

for(i in 1:nIter){
  setTxtProgressBar(pb, i)
  init =  sample(1:n)
  init = c(init, init[1])
  res <- optim(init, obj.fun.phase2, swap, method="SANN",
    control =
      list(maxit = 100000, temp = y2 , tmax = 1000, trace = TRUE,
                                                      REPORT = 10, fnscale = -1))

  newind = res$par[-length(res$par)]

  #Save all the 100 indexes for comparing later
  indList = rbind(indList, newind)
}
close(pb)
indList = indList[-1,]

#compare the average efficiency factor of animal to run*tag
z1.aveEff = apply(indList, 1, function(x) test(X.trt = Z1.mat[x,], blk.proj, Rep=Z1.rep)$ave.eff)

z1.aveEff[which(z1.aveEff == max(z1.aveEff))]

z1.canEff = apply(indList, 1, function(x) table(test(X.trt = Z1.mat[x,],  (mI(n) - Pb), Rep=Z1.rep)$can.eff))

z1.sd.canEff = apply(indList, 1, function(x) sd(test(X.trt = Z1.mat[x,], blk.proj, Rep=Z1.rep)$can.eff))

which(z1.aveEff == max(z1.aveEff[-which(z1.aveEff == max(z1.aveEff))]))

which(z1.aveEff == max(z1.aveEff[-c(which(z1.aveEff == max(z1.aveEff[-which(z1.aveEff == max(z1.aveEff))])),
which(z1.aveEff == max(z1.aveEff)))]))

PP = blk.proj %*% Z1.mat[bestInd,] %*% (mI(nZ1) - mK(nZ1)) %*%
     invInfMat((mI(n) - Pb), Z1.mat[bestInd,], mI(nZ1) - mK(nZ1)) %*%
      (mI(nZ1) - mK(nZ1)) %*% t(Z1.mat[bestInd,]) %*% blk.proj

test(X.trt = X.trt[bestInd,], blk.proj = PP, Rep=trt.rep )

new.Z1.mat=  Z1.mat[newInd,]

##############################################################################
#  Write a function that compute the remaining efficiency factor after eliminating the confounding
##############################################################################
trt.confounding = function(x){
  #x = indList[3,]
  new.Z1.mat =  Z1.mat[x,]
  colnames(new.Z1.mat) = sort(levels(interaction(LETTERS[1:nZ1])))
  new.Z1.des = apply(new.Z1.mat, 1,  function(x)
  colnames(new.Z1.mat)[which(as.logical(x))])
   new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))
  
  design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))
  
  design.df = cbind(design.df,
              Ani = t(new.Z1.des ),
              Trt = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                            as.character(phase1DesignEX1$Ani)),]$Trt)
  summary.aov.twoPhase(design.df,  blk.str2 = "Run",  
                        blk.str1 = "Ani", trt.str = "Tag + Trt")
                        
  trt.EF = as.character(summary.aov.twoPhase(design.df,  blk.str2 = "Run",  
                        blk.str1 = "Ani", trt.str = "Tag + Trt")$EF[,"eff.Trt"])
  
  numeric.trt.EF = sapply(strsplit(trt.EF, "/"), function(x) as.numeric(x)[1]/as.numeric(x)[2])
  
  numeric.trt.EF = numeric.trt.EF[-which(is.na(numeric.trt.EF))]
  
  #if((length(numeric.trt.EF) == 0) || (length(numeric.trt.EF) == 2)){
  #  return(0)
  #} else {
    return(numeric.trt.EF)
  #}
}

apply(indList, 1, trt.confounding)

trt.aveEff =  apply(indList, 1, trt.confounding)

trt.aveEff[which(trt.aveEff == max(trt.aveEff))]

trt.aveEff[which(trt.aveEff == max(trt.aveEff[-which(trt.aveEff == max(trt.aveEff))]))]

#compare the average efficiency factor of trt to (run*tag)/ani
trt.aveEff = apply(indList, 1, function(x) test(X.trt = X.trt[x,],
          blk.proj = blk.proj %*% Z1.mat[x,] %*% (mI(nZ1) - mK(nZ1)) %*%
     invInfMat( blk.proj, Z1.mat[x,], mI(nZ1) - mK(nZ1)) %*%
      (mI(nZ1) - mK(nZ1)) %*% t(Z1.mat[x,]) %*% blk.proj,
      Rep=trt.rep)$ave.eff)

trt.aveEff = apply(indList, 1, function(x) test(X.trt = X.trt[x,],
          blk.proj = (mI(n) - Pb) %*% Z1.mat[x,] %*% (mI(nZ1) - mK(nZ1)) %*%
     invInfMat( (mI(n) - Pb) , Z1.mat[x,], mI(nZ1) - mK(nZ1)) %*%
      (mI(nZ1) - mK(nZ1)) %*% t(Z1.mat[x,]) %*%(mI(n) - Pb) ,
      Rep=trt.rep)$ave.eff)

tag.aveEff = apply(indList, 1, function(x) test(X.trt = Zt,
          blk.proj = (mI(n) - Pb) %*% Z1.mat[x,] %*% (mI(nZ1) - mK(nZ1)) %*%
     invInfMat( (mI(n) - Pb) , Z1.mat[x,], mI(nZ1) - mK(nZ1)) %*%
      (mI(nZ1) - mK(nZ1)) %*% t(Z1.mat[x,]) %*%(mI(n) - Pb) ,
      Rep=trt.rep)$ave.eff)

 which(tag.aveEff == min(tag.aveEff))

trt.aveEff[which(trt.aveEff == max(trt.aveEff))]

test(X.trt = X.trt[indList[which(trt.aveEff == max(trt.aveEff)),],], blk.proj = PP, Rep=trt.rep )

test(X.trt = Z1.mat[indList[which(trt.aveEff == max(trt.aveEff)),],],
      blk.proj, Rep=Z1.rep)$ave.eff

new.Z1.mat =  Z1.mat[indList[which(trt.aveEff == max(trt.aveEff)),],]

new.Z1.mat =  Z1.mat[indList[5,],]


test(X.trt = new.Z1.mat,  (mI(n) - Pb), Rep=Z1.rep )
test(X.trt = new.Z1.mat,  blk.proj, Rep=Z1.rep )

#Confirm it with the cpp function I wrote
xxxx = cpp.new.fact.s.optimised(Z1.mat, blk.proj, nIter = 100)

test(X.trt = xxxx$b,  (mI(n) - Pb), Rep=Z1.rep)


test(X.trt = xxxx$b, blk.proj, Rep=Z1.rep)$ave.eff
test(X.trt = xxxx$b, blk.proj, Rep=Z1.rep)$can.eff

test(X.trt = xxxx$b, (mI(n) - Pb), Rep=Z1.rep)

#new.Z1.mat = xxxx$b


colnames(X.trt) = sort(levels(interaction(LETTERS[1:nTrt])))

trt.des = phase1DesignEX1[match(as.character(t(Z1.des )),
                          as.character(phase1DesignEX1$Ani)),]$Trt
X.trt = matrix(0, ncol = nTrt, nrow= nTrt*trt.rep )
X.trt[cbind(1:( nTrt*trt.rep ), trt.des)] = 1

new.Z1.mat =  Z1.mat[indList[2,],]

colnames(new.Z1.mat) = sort(levels(interaction(LETTERS[1:nZ1])))
new.Z1.des = apply(new.Z1.mat, 1,  function(x)
colnames(new.Z1.mat)[which(as.logical(x))])
 new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))

design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

design.df = cbind(design.df,
            Ani = t(new.Z1.des ),
            Trt = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                          as.character(phase1DesignEX1$Ani)),]$Trt)


matrix(design.df$Ani, nrow = nBlk, ncol = nPlot, byrow = TRUE)
(N = with(design.df, table(Ani, Run)))

N %*% t(N)

(N = with(design.df, table(Ani, Tag)))
N %*% t(N)

(N = with(design.df, table(Trt, Run)))
N %*% t(N)

(N = with(design.df, table(Trt, Tag)))
N %*% t(N)

trt.mat = matrix(design.df$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE)
trt.mat[which(trt.mat == "a")] = "controlled"
trt.mat[which(trt.mat == "b")] = "diseased"
trt.mat[which(trt.mat == "c")] = "treated"
trt.mat

summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", 
trt.str = "Tag + Trt")

summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", 
trt.str = "Trt + Tag")

summary.aov.twoPhase(design.df,  blk.str2 = "Run + Tag", blk.str1 = "Ani", 
trt.str = "Trt")

summary.aov.twoPhase(design.df,  blk.str2 = "Run+ Tag", blk.str1 = "Trt", 
trt.str = "Trt")

c(nTrt, nAni, nPlot, n/nPlot)

fileName = paste("design_", nTrt, "trt_", nAni, "ani_", nPlot, "tag_", n/nPlot, 
"run.Rdata", sep = "")

fileName

save(design.df, file = fileName)






design.summary = function(design.df, trtFirst = TRUE){

  n =  nrow(design.df)
  nBlk = nlevels(design.df$Run)
  nPlot = nlevels(design.df$Tag)
  Run.mat = with(design.df, as.matrix(table(1:n, Run)))
  Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
  Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
  Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))

  Pb = projMat(Run.mat)
  Pb1 = projMat(Tag.mat)

  blk.proj = (mI(n) - Pb) %*% (mI(n) - Pb1)

  print("Animal design:")
  print(matrix(design.df$Ani, nrow = nBlk, ncol = nPlot, byrow = TRUE))

  print("Animal incidence matrix:")
  print((N = with(design.df, table(Ani, Run))))

  print("Animal concurrence matrix:")
  print(N %*% t(N))

  print(test(X.trt = Ani.mat,  blk.proj, Rep=n/ncol(Ani.mat)))

  print("Treatment design:")
  print(matrix(design.df$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE))

  print("Treatment incidence matrix:")
  print((N = with(design.df, table(Trt, Run))))

  print("Treatment concurrence matrix:")
  print(N %*% t(N))

  PP = blk.proj %*% Ani.mat %*% (mI(ncol(Ani.mat)) - mK(ncol(Ani.mat))) %*%
       invInfMat(blk.proj, Ani.mat, mI(ncol(Ani.mat)) - mK(ncol(Ani.mat))) %*%
        (mI(ncol(Ani.mat)) - mK(ncol(Ani.mat))) %*% t(Ani.mat) %*% blk.proj

  print(test(X.trt = Trt.mat,  blk.proj = PP, Rep=n/ncol(Trt.mat)))

  if(trtFirst){
    summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Trt + Tag")
  }else{
    summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")

  }
}

design.summary(design.df)
