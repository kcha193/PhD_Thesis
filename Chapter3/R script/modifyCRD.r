sourceDir <- function(path, trace = TRUE, ...) {
  library(MASS)
  library(lattice)
  library(Rcpp)
  library(inline)
  library(RcppArmadillo)
  library(compiler)
  library(formatR)
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
  }
}



test =
function(X.trt, blk.proj, C.trt.mat = diag(ncol(X.trt)), Rep){
  if( length(dim(X.trt))==3) X.trt = X.trt[,,1]

  info.mat = C.trt.mat %*% t(X.trt) %*% blk.proj %*% X.trt %*% C.trt.mat
  trace = tr(info.mat)
  e.va = eigen(info.mat)$va
  e.vec =  eigen(info.mat)$vec
  can.eff = e.va[-which(e.va<1e-7)]/Rep
  list( trace = trace,
        nCan = length(can.eff),
        can.eff = can.eff,
        ave.eff =  1/mean(1/can.eff),
        e.vec = e.vec )
}

sourceDir(path = "C:/Users/kcha193/My Dropbox/R functions for two-phase experiments")

################################################################################

nTrt = 2
bRep = 9
tRep = 4
n = nTrt * bRep * tRep
nPlot = 8
nBlk = n/nPlot

nAni = nTrt * bRep

c(nTrt, nAni, nPlot, n/nPlot)

phase1DesignEX1 <- local({
    Ani = 1:nAni
    Trt = 1:nTrt
   
    data.frame(cbind(Ani, Trt))
})
phase1DesignEX1$Ani = as.factor(LETTERS[phase1DesignEX1$Ani])
phase1DesignEX1$Trt = as.factor(letters[phase1DesignEX1$Trt])

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



# Parameter's of treatment
trt.rep = n/nTrt

trt.des = as.numeric(phase1DesignEX1$Trt[match(LETTERS[Z1.des], phase1DesignEX1$Ani)])
X.trt = matrix(0, ncol = nTrt, nrow = nTrt * trt.rep)
X.trt[cbind(1:(nTrt * trt.rep), trt.des)] = 1 

################################################################################
swap.trt <- function(ind) {
    
    idx <- seq(2, length(ind) - 2)
    changepoints <- sample(idx, size = 2, replace = FALSE)
    #checks for omitting the swap of two indtical animals/ trtments
    check1 = all(X.trt [ind[changepoints[1]],] == X.trt[ind[changepoints[2]],]) 
    
    while(check1){
          changepoints <- sample(idx, size = 2, replace = FALSE)
          check1 = all(X.trt [ind[changepoints[1]],] == X.trt [ind[changepoints[2]],])   
    }
    tmp <- ind[changepoints[1]]
    ind[changepoints[1]] <- ind[changepoints[2]]
    ind[changepoints[2]] <- tmp
    return(ind)
}


#Objective functions 
obj.fun.trt = function(ind) {
         
    Rep = trt.rep
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*% blk.proj %*% newX
    e.va = eigen(info.mat, only.values = TRUE)$va
    can.eff = e.va[-which(e.va < 1e-07)]/Rep

    aveEff = 1/mean(1/can.eff)

    return(aveEff)
}


#Objective functions 
obj.fun.z1 = function(ind) {
     
    Rep = Z1.rep
    newZ1 = Z1.mat[ind[-length(ind)], ]
    info.mat = t(newZ1) %*% blk.proj %*% newZ1
    e.va = eigen(info.mat, only.values = TRUE)$va
    can.eff = e.va[-which(e.va < 1e-07)]/Rep
    
    aveEff = 1/mean(1/can.eff)
       
    return(aveEff)
}

################################################################################
#Assigning treatments to Phase 2 block structure

init = sample(1:n)
newInit = c(init, init[1])

y2 = 1

res <- optim(newInit, obj.fun.trt, swap.trt, method="SANN",
      control =
          list(maxit = 50000  , temp = y2 , tmax = 1000, trace = TRUE, 
                                                  REPORT = 10, fnscale = -1))

newind = res$par[-length(res$par)]

trt.des[newind]


################################################################################
#Assigning Phase 1 block to Phase 2 block structure


Z1.des = as.factor(rep(1:nZ1, Z1.rep))

table.Z1.des = table(Z1.des)
Z1.des = numeric(n)
numeric.Z1.des = as.numeric(phase1DesignEX1$Ani)
trt.factor = phase1DesignEX1$Trt

for (i in 1:n) {
    
    Z1.des[i] = sample(names(table.Z1.des)[which(letters[trt.des[newind][i]] == trt.factor)])[1]
    
    table.Z1.des[Z1.des[i]] = table.Z1.des[Z1.des[i]] - 1
    
    if (any(table.Z1.des == 0)) {
        trt.factor = trt.factor[-as.numeric(which(table.Z1.des == 0))]
        table.Z1.des = table.Z1.des[-which(table.Z1.des == 0)]
        
    }
    
} 

Z1.des = as.numeric(Z1.des)
Z1.mat = matrix(0, ncol = nZ1, nrow = nZ1 * Z1.rep)
Z1.mat[cbind(1:(nZ1 * Z1.rep), Z1.des)] = 1




swap.z1 <- function(ind) {
    
    idx <- seq(2, length(ind) - 2)
    changepoints <- sample(idx, size = 2, replace = FALSE)
    #checks for omitting the swap of two indtical animals/ trtments
    check1 = all(Z1.mat[ind[changepoints[1]],] == Z1.mat[ind[changepoints[2]],]) 
    
    check2 =  phase1DesignEX1[which(Z1.mat[ind[changepoints[1]],]==1),]$Trt ==  
                    phase1DesignEX1[which(Z1.mat[ind[changepoints[2]],]==1),]$Trt
    
    while(check1 || check2){
          changepoints <- sample(idx, size = 2, replace = FALSE)
          check1 = all(Z1.mat[ind[changepoints[1]],] == Z1.mat[ind[changepoints[2]],]) 
          check2 =  phase1DesignEX1[which(Z1.mat[ind[changepoints[1]],]==1),]$Trt ==  
                    phase1DesignEX1[which(Z1.mat[ind[changepoints[2]],]==1),]$Trt
  
    }
    tmp <- ind[changepoints[1]]
    ind[changepoints[1]] <- ind[changepoints[2]]
    ind[changepoints[2]] <- tmp
    return(ind)
}



#Objective functions 
obj.fun.z1 = function(ind) {
     
    Rep = Z1.rep
    newZ1 = Z1.mat[ind[-length(ind)], ]
    info.mat = t(newZ1) %*% blk.proj %*% newZ1
    e.va = eigen(info.mat, only.values = TRUE)$va
    can.eff = e.va[-which(e.va < 1e-07)]/Rep
    
    aveEff = 1/mean(1/can.eff)
       
    return(aveEff)
}



init = 1:n #sample(1:n)
newInit = c(init, init[1])

y2 = 0.2

res <- optim(newInit, obj.fun.z1, swap.z1, method="SANN",
      control =
          list(maxit = 100000  , temp = y2 , tmax = 1000, trace = TRUE, 
                                                  REPORT = 10, fnscale = -1))

newind = res$par[-length(res$par)]

test(X.trt = Z1.mat[newind,], mI(n) - Pb, Rep=Z1.rep)$ave.eff
test(X.trt = Z1.mat[newind,], mI(n) - Pb, Rep=Z1.rep)$can.eff

new.Z1.mat=  Z1.mat[newind,]


################################################################################
#Set the design from the search 

 colnames(new.Z1.mat) = sort(levels(interaction(LETTERS[1:nZ1])))
new.Z1.des = apply(new.Z1.mat, 1,  function(x)
colnames(new.Z1.mat)[which(as.logical(x))])
 new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))

design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

design.df = cbind(design.df,
            Ani = t(new.Z1.des ),
            Trt = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                          as.character(phase1DesignEX1$Ani)),]$Trt)

################################################################################
#Theortical ANOVA table 

summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", 
trt.str = "Tag + Trt")
