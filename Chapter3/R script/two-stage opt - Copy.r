sourceDir <- function(path, trace = TRUE, ...) {
  library(MASS)
  library(lattice)
  library(inline)
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

multiLETTERS = function(x) {

    if (max(x) < 26) {
        LETTERS[x]
    } else {
        newLETTERS = sort(sapply(strsplit(levels(interaction(LETTERS, LETTERS)), "\\."), function(x) paste(x, collapse = "", sep = "")))
        c(LETTERS[1:26], newLETTERS[x[-(1:26)] - 26])
    }
}

multiLetters = function(x) {

    if (max(x) < 26) {
        letters[x]
    } else {
        newLETTERS = sort(levels(interaction(letters, letters)))
        c(letters[1:26], newLETTERS[x[-(1:26)] - 26])
    }


}


is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


newLETTERS = sort(sapply(strsplit(levels(interaction(LETTERS, LETTERS)), "\\."), function(x) paste(x, collapse = "", sep = "")))

sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R")
sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/optimTE/R")

################################################################################

nTrt = 4
bRep = 3
tRep = 6
n = nTrt * bRep * tRep
nPlot = 6
nBlk = n/nPlot
trt.rep = n/nTrt

nAni = nTrt * bRep

c(nTrt, nAni, nPlot, n/nPlot)

phase1DesignEX1 <- local({
    Ani = 1:nAni
    Pro = 1:nTrt
    Sam = rep(1:tRep, each = nAni)
    Inter = interaction(Pro, Sam)
    data.frame(cbind(Ani, Pro, Sam, Inter))
})
    phase1DesignEX1$Ani = as.factor(multiLetters(phase1DesignEX1$Ani))
    phase1DesignEX1$Pro = as.factor(multiLetters(phase1DesignEX1$Pro, FALSE))
    phase1DesignEX1$Sam = as.factor(phase1DesignEX1$Sam)
    phase1DesignEX1$Inter =  interaction(phase1DesignEX1$Pro, phase1DesignEX1$Sam, sep = "")
    
summaryAovOnePhase(phase1DesignEX1, blk.str = "Ani", trt.str = "Pro * Sam")


# Parameter's of block
nBlk = n/nPlot

Zb = matrix(0, ncol = nBlk, nrow = n)
Z.des = rep(1:nBlk, each = nPlot)
Zb[cbind(1:n, Z.des)] <- 1


Pb = projMat(Zb)

withBlock = (identityMat(n) - Pb) 

blk.proj = withBlock


Z1.des = initialCRD(nTrt = nTrt, bRep = bRep, tRep = tRep, nPlot = nPlot)

# Parameter's of block structure of Phase 1
   nZ1 = nrow(phase1DesignEX1)
    Z1.rep = n/nZ1

Z1.mat = matrix(0, ncol = nZ1, nrow = nZ1 * Z1.rep)
Z1.mat[cbind(1:(nZ1 * Z1.rep), Z1.des)] = 1

nTrt = nlevels(phase1DesignEX1$Inter)

# Parameter's of treatment
trt.rep = n/nTrt

trt.des = as.numeric(phase1DesignEX1$Inter[match(multiLetters(Z1.des), phase1DesignEX1$Ani)])
X.trt = matrix(0, ncol = nTrt, nrow = nTrt * trt.rep)
X.trt[cbind(1:(nTrt * trt.rep), trt.des)] = 1 

 
#################################################################################

getPair = function(x, n) {
    y = ceiling(x/nPlot)
    if(nPlot == 4){
    if ((n/nPlot)%%2 == 1) {
         
        if (x == n) {
            c(x, x - 1)
        } else if (x == (n - 1)) {
            c(x, x + 1)
        } else if (x == (n - 2)) {
            c(x, x - 1)
        } else if (x == (n - 3)) {
            c(x, x + 1)
        } else if (y%%2 == 0) {
            if (x%%2 == 1) {
                c(x, x - 3)
            } else {
                c(x, x - 5)
            }
        } else {
            if (x%%2 == 0) {
                c(x, x + 3)
            } else {
                c(x, x + 5)
            }
        }
    } else {
        if (y%%2 == 0) {
            if (x%%2 == 1) {
                c(x, x - 3)
            } else {
                c(x, x - 5)
            }
        } else {
            if (x%%2 == 0) {
                c(x, x + 3)
            } else {
                c(x, x + 5)
            }
        }
    }
    } else {
    if ((n/nPlot)%%2 == 1) {
         
        if (x == n) {
            c(x, x - 1)
        } else if (x == (n - 1)) {
            c(x, x + 1)
        } else if (x == (n - 2)) {
            c(x, x - 1)
        } else if (x == (n - 3)) {
            c(x, x + 1)
        } else if (x == (n - 4)) {
            c(x, x - 1)
        } else if (x == (n - 5)) {
            c(x, x + 1)
        } else if (x == (n - 6)) {
            c(x, x - 1)
        } else if (x == (n - 7)) {
            c(x, x + 1)
        } else if (y %% 2 == 0) {
            if (x%%2 == 1) {
                c(x, x - 7)
            } else {
                c(x, x - 9)
            }
        } else {
            if (x%%2 == 0) {
                c(x, x + 7)
            } else {
                c(x, x + 9)
            }
        }
    } else {
        if (y%%2 == 0) {
            if (x%%2 == 1) {
                c(x, x - 7)
            } else {
                c(x, x - 9)
            }
        } else {
            if (x%%2 == 0) {
                c(x, x + 7)
            } else {
                c(x, x + 9)
            }
        }
    }
    
    
    }
} 
 

swap.new <- function(ind) {
    
    idx <- seq(2, length(ind) - 2)
    changepoints <- sample(idx, size = 2, replace = FALSE)
    #Two checks for omitting the swap of two indtical animals and animals in the same block
    check1 = all(Z1.mat[ind[changepoints[1]],] == Z1.mat[ind[changepoints[2]],]) 
    
    while(check1){
          changepoints <- sample(idx, size = 2, replace = FALSE)
          check1 = all(Z1.mat[ind[changepoints[1]],] == Z1.mat[ind[changepoints[2]],])   
    }
    tmp <- ind[getPair(changepoints[1], length(ind) - 1)]
    ind[getPair(changepoints[1], length(ind) - 1)] <- ind[getPair(changepoints[2], length(ind) - 1)]
    ind[getPair(changepoints[2], length(ind) - 1)] <- tmp
    return(ind)
}

swap.stage1.new <- function(ind) {
    
    idx <- seq(2, length(ind) - 2)
    changepoints <- sample(idx, size = 2, replace = FALSE)
    #Two checks for omitting the swap of two indtical animals and animals in the same block
    check1 = all(Z1.mat[ind[changepoints[1]],] == Z1.mat[ind[changepoints[2]],]) 
    check2 = diff((changepoints - 1)%/%nPlot) != 0
    
    while(check1 || check2){
          changepoints <- sample(idx, size = 2, replace = FALSE)
          check1 = all(Z1.mat[ind[changepoints[1]],] == Z1.mat[ind[changepoints[2]],]) 
          check2 = diff((changepoints - 1)%/%nPlot) != 0

    }
       
    tmp <- ind[getPair(changepoints[1], length(ind) - 1)]
    ind[getPair(changepoints[1], length(ind) - 1)] <- ind[getPair(changepoints[2], length(ind) - 1)]
    ind[getPair(changepoints[2], length(ind) - 1)] <- tmp
    return(ind)
}

swap.stage2.new <- function(ind) {
    
    idx <- seq(2, length(ind) - 2)
    changepoints <- sample(idx, size = 2, replace = FALSE)
    #Two checks for omitting the swap of two indtical animals and animals in the same block
    check1 = all(Z1.mat[ind[changepoints[1]],] == Z1.mat[ind[changepoints[2]],]) 
    check2 = diff(changepoints%%nPlot) != 0
    
    while(check1 || check2){
          changepoints <- sample(idx, size = 2, replace = FALSE)
          check1 = all(Z1.mat[ind[changepoints[1]],] == Z1.mat[ind[changepoints[2]],]) 
          check2 = diff(changepoints%%nPlot) != 0
     }
    tmp <- ind[getPair(changepoints[1], length(ind) - 1)]
    ind[getPair(changepoints[1], length(ind) - 1)] <- ind[getPair(changepoints[2], length(ind) - 1)]
    ind[getPair(changepoints[2], length(ind) - 1)] <- tmp
    return(ind)
}


swap.stage3.new <- function(ind) {

     idx <- seq(2, length(ind) - 2)
    changepoints <- sample(idx, size = 2, replace = FALSE)
    #Two checks for omitting the swap of two indtical animals and animals in the same block
    check1 = all(Z1.mat[ind[changepoints[1]],] == Z1.mat[ind[changepoints[2]],])
    check2 = diff((changepoints - 1)%/%nPlot) == 0
    check3 = diff(changepoints%%nPlot) == 0

    while(check1 || check2 || check3){
          changepoints <- sample(idx, size = 2, replace = FALSE)
          check1 = all(Z1.mat[ind[changepoints[1]],] == Z1.mat[ind[changepoints[2]],])
          check2 = diff((changepoints - 1)%/%nPlot) == 0
          check3 = diff(changepoints%%nPlot) == 0
    }
    tmp <- ind[getPair(changepoints[1], length(ind) - 1)]
    ind[getPair(changepoints[1], length(ind) - 1)] <- ind[getPair(changepoints[2], length(ind) - 1)]
    ind[getPair(changepoints[2], length(ind) - 1)] <- tmp
    return(ind)
}

 
C.pro  = (mI(nCag) - mK(nCag)) %x%  mK(nAni/nCag) 
C.ani =   mI(nCag)%x%  (mI(nAni/nCag) - mK(nAni/nCag))
cage.Rep = n/nCag

C.ani =   mI(nCag)%x%  (mI(nAni/nCag) - mK(nAni/nCag))
ani.Rep = n/nAni
trt.Rep = trt.rep

   

#Objective functions 
try.obj.fun.CRD <-
function(ind) {
        
        Rep = Z1.rep
        newZ1 = Z1.mat[ind[-length(ind)], ]
        info.mat = t(newZ1) %*% blk.proj %*% newZ1
        e.va = eigen(info.mat, only.values = TRUE)$va
        can.eff = e.va[-which(e.va < 1e-07)]/Rep
        
        aveEff = 1/mean(1/can.eff)
        
        Rep = trt.rep
        newX = X.trt[ind[-length(ind)], ]
        info.mat = t(newX) %*% blk.proj %*% newX
        e.va = eigen(info.mat, only.values = TRUE)$va
        can.eff = e.va[-which(e.va < 1e-07)]/Rep
        len = length(can.eff)
        
        if(len == 0) return(0)
        
        aveEff = (3/4) * aveEff + 1/4 * ((len + 1/mean(1/can.eff))/nTrt)
        
        return(aveEff)
    }


obj.fun.CRD = function(ind) {

  ans = try(try.obj.fun.CRD(ind), silent = TRUE)
      
  ifelse(class(ans) == "try-error", 0, ans)
}



################################################################################
#Compare different optimality criteria

set.seed(527)
init = sample(1:n)

newInit = c(init, init[1])
old.time = proc.time()

res <- optim(newInit, obj.fun.new, swap, method = "SANN", control = list(maxit = 1e+6, 
 trace = TRUE, REPORT = 1000, fnscale = -1))

proc.time() - old.time

old.time = proc.time()

set.seed(820)
init = sample(1:n)

newInit = c(init, init[1])

res <- optim(newInit, obj.fun.ms, swap, method = "SANN", control = list(maxit = 1e5, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

proc.time() - old.time


res <- optim(newInit, obj.fun.s, swap, method = "SANN", control = list(maxit = 1e5, 
    temp = 10, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

old.time = proc.time()

set.seed(527)
init = sample(1:n)

newInit = c(init, init[1])

res <- optim(newInit, obj.fun.ave, swap, method = "SANN", control = list(maxit = 1e5, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

proc.time() - old.time

obj.fun.ms(newInit)

old.time = proc.time()

set.seed(527)
init = sample(1:n)

newInit = c(init, init[1])

res <- optim(newInit, obj.fun.old, swap, method = "SANN", control = list(maxit = 1e+5, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

proc.time() - old.time


set.seed(527)
init = sample(1:n)

newInit = c(init, init[1])

res <- optim(newInit, obj.fun.old1, swap, method = "SANN", control = list(maxit = 1e+5, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

old.time = proc.time()


################################################################################
#Compare random starting temperature versus modified starting temperature
set.seed(1001)
init = sample(1:n)

newInit = c(init, init[1])

obj.fun.CRD(newInit) 

old.time = proc.time()
res <- optim(newInit, obj.fun.CRD, swap, method = "SANN", control = list(maxit = 1e+5, 
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
proc.time() - old.time


initialTemp = function(iter, newInit, swap, obj.fun){
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  U = numeric(iter)
  temp = matrix(0, nrow = iter, ncol = length(newInit))
  temp[1,] = newInit
  U[1] =   obj.fun(temp[1,])
  tol = 10
  check.U = numeric(tol) 
  
  for(i in 2:iter){
    setTxtProgressBar(pb, i)
    temp[i,] = swap(temp[i-1,])
    U[i] = obj.fun(temp[i,])
    
    check.U = c(check.U, U[i])
    check.U = check.U[-1]
    
    if(all(outer(check.U,check.U,"=="))) temp[i,] = newInit
  }
  close(pb)
  
  if(all(U==0) || all(diff(U)<1e-8)){
    return(list(y2 = 100, y2low = 1e-5, newInit = newInit))
  
  }else{
    xxx = max(U, na.rm = TRUE) - range(U[which(U < (max(U, na.rm = TRUE)-1e-10))], na.rm = TRUE)
    
    if(isTRUE(all.equal(xxx[1], xxx[2]))) xxx[2] = 1e-5
     
    newInit = temp[which(U == max(U))[1],]
    return(list(y2 = max(xxx), y2low = min(xxx), newInit = newInit))
  }
}   

newInit = c(init, init[1])
obj.fun.CRD(newInit) 

temp = initialTemp(1e4, newInit, swap = swap, obj.fun = obj.fun.CRD)
temp 
(y2 = temp$y2)

newInit = temp$newInit 

res <- optim(newInit, obj.fun.CRD, swap, method = "SANN", control = list(maxit = 9e4, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

temp$y2/ log(((t-1) %/% 1000)*1000 + exp(1))

exp(temp.save$y2/temp.save$y2low)

10/ log(((100000-1) %/% 1000)*1000 + exp(1))

(10/7/7/7/7/7/7/7/7/7)/log(((10000-1) %/% 1000)*1000 + exp(1))

################################################################################
#Compare	The accelerated cooling method versus standard cooling method.

set.seed(1001)
init = sample(1:n)
newInit = c(init, init[1])
obj.fun.CRD(newInit)
  
y2 = 10

res <- optim(newInit, obj.fun.CRD, swap, method = "SANN", control = list(maxit = 1e+4, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

for (i in 1:9) {
    y2 = y2/7
    print(y2)
    res <- optim(res$par, obj.fun.CRD, swap, method = "SANN", control = list(maxit = 1e+4, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

}


set.seed(1001)
init = sample(1:n)
newInit = c(init, init[1])

temp =  initialTemp(iter = 1e4, newInit = newInit, swap = swap, obj.fun = obj.fun.CRD) 
temp
(y2 = temp$y2)

newInit = temp$newInit 

inter = exp(log(y2/temp$y2low)/8)

res <- optim(newInit, obj.fun.CRD, swap, method = "SANN", control = list(maxit = 1e+4, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

for (i in 1:8) {
    y2 = y2/inter
    print(y2)
    
    res <- optim(res$par, obj.fun.CRD, swap, method = "SANN", control = list(maxit = 1e+4, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
}


################################################################################
#Modified accelerated cooling methods

set.seed(1001)
init = sample(1:n)
newInit = c(init, init[1])
iter = 1e4
  
y2 = 10

old = obj.fun.CRD(newInit)

res <- optim(newInit, obj.fun.CRD, swap, method = "SANN", control = list(maxit = iter, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

cur = obj.fun.CRD(res$par)

while ((cur - old) > 1e-06) {
  old = cur
  
  res <- optim(res$par, obj.fun.CRD, swap, method = "SANN", control = list(maxit =iter, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

  cur = obj.fun.CRD(res$par)

}


for (i in 1:9) {
    y2 = y2/7
    print(y2)
    
    old = obj.fun.CRD(res$par)
    res <- optim(res$par, obj.fun.CRD, swap, method = "SANN", control = list(maxit = iter, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
    
    cur = obj.fun.CRD(res$par)
  
  while ((cur - old) > 1e-06) {
    old = cur
    
    res <- optim(res$par, obj.fun.CRD, swap, method = "SANN", control = list(maxit = iter, 
      temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
  
    cur = obj.fun.CRD(res$par)
  
  }
}

##############################################################################

set.seed(1001)
init = sample(1:n)
newInit = c(init, init[1])
iter = 1e4

(old = obj.fun.CRD(newInit))

temp =  initialTemp(iter = iter, newInit = newInit, swap = swap, obj.fun = obj.fun.CRD) 

newInit = temp$newInit 

cur = obj.fun.CRD(temp$newInit)

while ((cur - old) > 1e-06) {
  old = cur
  
  temp =  initialTemp(iter = iter, newInit = newInit, swap = swap, obj.fun = obj.fun.CRD) 

  cur = obj.fun.CRD(temp$newInit)

}

(y2 = temp$y2)

newInit = temp$newInit 

old = obj.fun.CRD(newInit)

res <- optim(newInit, obj.fun.CRD, swap, method = "SANN", control = list(maxit = iter, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

cur = obj.fun.CRD(res$par)

while ((cur - old) > 1e-06) {
  old = cur
  
  res <- optim(res$par, obj.fun.CRD, swap, method = "SANN", control = list(maxit =iter, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

  cur = obj.fun.CRD(res$par)

}

inter = exp(log(y2/temp$y2low)/8)

for (i in 1:8) {
    y2 = y2/inter
    print(y2)
    
    old = obj.fun.CRD(res$par)
    res <- optim(res$par, obj.fun.CRD, swap, method = "SANN", control = list(maxit = iter, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
    
    cur = obj.fun.CRD(res$par)
  
  while ((cur - old) > 1e-06) {
    old = cur
    
    res <- optim(res$par, obj.fun.CRD, swap, method = "SANN", control = list(maxit = iter, 
      temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
  
    cur = obj.fun.CRD(res$par)
  
  }    

}

################################################################################
#Compare random starting design versus modified starting design

init = 1:n #sample(1:n)
newInit = c(init, init[1])
iter = 1e4

(old = obj.fun.CRD(newInit))

temp =  initialTemp(iter = iter, newInit = newInit, swap = swap, obj.fun = obj.fun.CRD) 

newInit = temp$newInit 

cur = obj.fun.CRD(temp$newInit)

while ((cur - old) > 1e-06) {
  old = cur
  
  temp =  initialTemp(iter = iter, newInit = newInit, swap = swap, obj.fun = obj.fun.CRD) 

  cur = obj.fun.CRD(temp$newInit)

}

(y2 = temp$y2)

newInit = temp$newInit 

old = obj.fun.CRD(newInit)

res <- optim(newInit, obj.fun.CRD, swap, method = "SANN", control = list(maxit = iter, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

cur = obj.fun.CRD(res$par)

while ((cur - old) > 1e-06) {
  old = cur
  
  res <- optim(res$par, obj.fun.CRD, swap, method = "SANN", control = list(maxit =iter, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

  cur = obj.fun.CRD(res$par)

}

inter = exp(log(y2/temp$y2low)/8)

for (i in 1:8) {
    y2 = y2/inter
    print(y2)
    
    old = obj.fun.CRD(res$par)
    res <- optim(res$par, obj.fun.CRD, swap, method = "SANN", control = list(maxit = iter, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
    
    cur = obj.fun.new(res$par)
  
  while ((cur - old) > 1e-06) {
    old = cur
    
    res <- optim(res$par, obj.fun.CRD, swap, method = "SANN", control = list(maxit = iter, 
      temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
  
    cur = obj.fun.CRD(res$par)
  
  }    

}

################################################################################
#Compare	Pair swapping method versus one-to-one swapping method.

init = 1:n #sample(1:n)
newInit = c(init, init[1])
iter = 1e4

(old = obj.fun.CRD(newInit))

temp =  initialTemp(iter = iter, newInit = newInit, swap = swap.new, obj.fun = obj.fun.CRD) 

newInit = temp$newInit 

cur = obj.fun.CRD(temp$newInit)

while ((cur - old) > 1e-06) {
  old = cur
  
  temp =  initialTemp(iter = iter, newInit = newInit, swap = swap.new, obj.fun = obj.fun.CRD) 

  cur = obj.fun.CRD(temp$newInit)

}

(y2 = temp$y2)

newInit = temp$newInit 

old = obj.fun.CRD(newInit)

res <- optim(newInit, obj.fun.CRD, swap.new, method = "SANN", control = list(maxit = iter, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

cur = obj.fun.CRD(res$par)

while ((cur - old) > 1e-06) {
  old = cur
  
  res <- optim(res$par, obj.fun.CRD, swap.new, method = "SANN", control = list(maxit =iter, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

  cur = obj.fun.CRD(res$par)

}

inter = exp(log(y2/temp$y2low)/8)

for (i in 1:8) {
    y2 = y2/inter
    print(y2)
    
    old = obj.fun.CRD(res$par)
    res <- optim(res$par, obj.fun.CRD, swap.new, method = "SANN", control = list(maxit = iter, 
    temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
    
    cur = obj.fun.CRD(res$par)
  
  while ((cur - old) > 1e-06) {
    old = cur
    
    res <- optim(res$par, obj.fun.CRD, swap.new, method = "SANN", control = list(maxit = iter, 
      temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
  
    cur = obj.fun.CRD(res$par)
  
  }    

}

################################################################################
#Compare Three-stages swapping method versus standard swapping method. 
set.seed(527)
init = sample(1:n)
newInit = c(init, init[1])

temp = initialTemp(11111, newInit, swap.stage1, obj.fun.new) 
(y2.stage1 = temp$y2)

res <- optim(temp$newInit, obj.fun.new, swap.stage1, method = "SANN", control = list(maxit = 3e5, 
    temp = y2.stage1, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

temp = initialTemp(11111, newInit, swap.stage2, obj.fun.new) 
(y2.stage2 = temp$y2)

init = 
if(obj.fun.new(temp$newInit)> obj.fun.new(res$par)){
  temp$newInit
} else {
  res$par 
}

res1 <- optim(init, obj.fun.new, swap.stage2, method = "SANN", control = list(maxit = 3e5,
    temp = y2.stage2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))


temp = initialTemp(11111, newInit, swap.stage3, obj.fun.new) 
(y2.stage3 = temp$y2)

init = 
if(obj.fun.new(temp$newInit)> obj.fun.new(res1$par)){
  temp$newInit
} else {
  res1$par 
}

res1 <- optim(temp$newInit, obj.fun.new, swap.stage3, method = "SANN", control = list(maxit = 3e5,
    temp = y2.stage3, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))


############################################################################################

test.obj.fun.CRD = function(ind) {

        
    #inforamtion in within runs and tags stratum while the cage information is eliminated
   
       Rep = Z1.rep
       newZ1 = Z1.mat[ind[-length(ind)], ]
       PP = projMat(newZ1) %*% blk.proj  
 
        Rep = trt.rep
       newX = X.trt[ind[-length(ind)], ]
       info.mat = t(newX) %*% blk.proj %*% newX
       e.va = eigen(info.mat, only.values = TRUE)$va
       can.eff = e.va[-which(e.va < 1e-07)]/Rep
       len = length(can.eff)

    
    return(list(aveEff = 100/mean(1/can.eff), res = tr(PP) - len))
}

############################################################################################################
#Final function

twoStageCRD =  function(init, iter = 10000, obj.fun, test.obj.fun.new) {
    print(c(nTrt, bRep, tRep, nPlot))
   
    newInit = c(init, init[1])

    old = obj.fun(newInit)

    cat("Level: 1, Finding the temperatures for both stages\n")

    temp = initialTemp(iter, newInit, swap.stage2, obj.fun)

    res.stage2 = temp$newInit
    y2.stage2 = temp$y2
    y2low.stage2 = temp$y2low


    temp = initialTemp(iter, newInit, swap.stage1, obj.fun)

    res.stage1 = temp$newInit
    y2.stage1 = temp$y2
    y2low.stage1 = temp$y2low

    
    cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
    cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
  
    test.newInit = apply(rbind(res.stage1,res.stage2),1, obj.fun)
    newInit = rbind(res.stage1,res.stage2)[which(test.newInit == max(test.newInit))[1],]
    
    cur = obj.fun(newInit)
    print(test.obj.fun.new(newInit))
    
    while ((cur - old) > 1e-06) {
        
        temp = initialTemp(iter, newInit, swap.stage2, obj.fun)
    
        res.stage2 = temp$newInit
        y2.stage2 = temp$y2
        y2low.stage2 = temp$y2low
    
    
        temp = initialTemp(iter, newInit, swap.stage1, obj.fun)
    
        res.stage1 = temp$newInit
        y2.stage1 = temp$y2
        y2low.stage1 = temp$y2low
    
        
           
        test.newInit = apply(rbind(res.stage1,res.stage2),1, obj.fun)
        newInit = rbind(res.stage1,res.stage2)[which(test.newInit == max(test.newInit))[1],]

        old = cur
        
        cur = obj.fun(newInit)
        print(test.obj.fun.new(newInit))     
   
        cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
        cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
     
     }
      

    old = obj.fun(newInit)
  
    cat("Level: 2, Current temp: ", y2.stage2, ",", y2.stage1, ",", y2, "\n")

    res <- optim(newInit, obj.fun, swap.stage2, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
        tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

    res1 <- optim(res$par, obj.fun, swap.stage1, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
        tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))


    cur = obj.fun(res1$par)

    print(test.obj.fun.new(res1$par))
    
    if(cur == 100) return(res1$par[-length(res1$par)])

    while ((cur - old) > 1e-06) {
        res1 <- optim(res1$par, obj.fun, swap.stage2, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        res1 <- optim(res1$par, obj.fun, swap.stage1, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        old = cur
        cur = obj.fun(res1$par)

        print(test.obj.fun.new(res1$par))


        if(cur == 100) return(res1$par[-length(res1$par)])
    }

    c(nTrt, bRep, tRep, nPlot)

    inter.stage1 = exp(log(y2.stage1/y2low.stage1)/8)
    inter.stage2 = exp(log(y2.stage2/y2low.stage2)/8)

    i = 1
    unstable = FALSE

    # for (i in 1:8) {
    while (i < 9 || unstable) {
        unstable = FALSE
        y2.stage2 = y2.stage2/inter.stage2
        y2.stage1 = y2.stage1/inter.stage1
        
        old = obj.fun(res1$par)

        cat("Level: ", i + 2, ", Current temp: ", y2.stage2, ",", y2.stage1, "\n")
        res1 <- optim(res1$par, obj.fun, swap.stage2, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        res1 <- optim(res1$par, obj.fun, swap.stage1, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
        
    
        cur = obj.fun(res1$par)
          print(test.obj.fun.new(res1$par))

        if(cur == 100) return(res1$par[-length(res1$par)])

        while ((cur - old) > 1e-06) {
            unstable = TRUE
            res1 <- optim(res1$par, obj.fun, swap.stage2, method = "SANN", control = list(maxit = iter,
                temp = y2.stage2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

            res1 <- optim(res1$par, obj.fun, swap.stage1, method = "SANN", control = list(maxit = iter,
                temp = y2.stage1, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

            old = cur
            cur = obj.fun(res1$par)
            print(test.obj.fun.new(res1$par))

            if(cur == 100) return(res1$par[-length(res1$par)])

        }


        i = i + 1
    }

    print(c(nTrt, bRep, tRep, nPlot))

    return(res1$par[-length(res1$par)])
}

set.seed(1001)
bestInd = twoStageCRD(init =  sample(1:n), iter = 5000, obj.fun = obj.fun.CRD, 
test.obj.fun.new = test.obj.fun.CRD)

############################################################################################################
#Final function

threeStageCRD =  function(init, iter = 10000, obj.fun, test.obj.fun.new) {
    print(c(nTrt, bRep, tRep, nPlot))
   
    newInit = c(init, init[1])

    old = obj.fun(newInit)

    cat("Level: 1, Finding the temperatures for both stages\n")

    temp = initialTemp(iter, newInit, swap.stage2, obj.fun)

    res.stage2 = temp$newInit
    y2.stage2 = temp$y2
    y2low.stage2 = temp$y2low


    temp = initialTemp(iter, newInit, swap.stage1, obj.fun)

    res.stage1 = temp$newInit
    y2.stage1 = temp$y2
    y2low.stage1 = temp$y2low


    temp = initialTemp(iter, newInit,swap.stage3, obj.fun)

    res = temp$newInit
    y2 = temp$y2
    y2low = temp$y2low
    
    cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
    cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
    cat("Stage 0 temperature range: ", y2, " to ", y2low, "\n")


    test.newInit = apply(rbind(res,res.stage1,res.stage2),1, obj.fun)
    newInit = rbind(res,res.stage1,res.stage2)[which(test.newInit == max(test.newInit))[1],]
    
    cur = obj.fun(newInit)
    print(test.obj.fun.new(newInit))
    
    while ((cur - old) > 1e-06) {
        
        temp = initialTemp(iter, newInit, swap.stage2, obj.fun)
    
        res.stage2 = temp$newInit
        y2.stage2 = temp$y2
        y2low.stage2 = temp$y2low
    
    
        temp = initialTemp(iter, newInit, swap.stage1, obj.fun)
    
        res.stage1 = temp$newInit
        y2.stage1 = temp$y2
        y2low.stage1 = temp$y2low
    
    
        temp = initialTemp(iter, newInit, swap.stage3, obj.fun)
    
        res = temp$newInit
        y2 = temp$y2
        y2low = temp$y2low
           
        test.newInit = apply(rbind(res,res.stage1,res.stage2),1, obj.fun)
        newInit = rbind(res,res.stage1,res.stage2)[which(test.newInit == max(test.newInit))[1],]

        old = cur
        
        cur = obj.fun(newInit)
        print(test.obj.fun.new(newInit))     
   
        cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
        cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
        cat("Stage 0 temperature range: ", y2, " to ", y2low, "\n")
    
     }
      

    old = obj.fun(newInit)
  
    cat("Level: 2, Current temp: ", y2.stage2, ",", y2.stage1, ",", y2, "\n")

    res <- optim(newInit, obj.fun, swap.stage2, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
        tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

    res1 <- optim(res$par, obj.fun, swap.stage1, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
        tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

    res1 <- optim(res1$par, obj.fun, swap.stage3, method = "SANN", control = list(maxit = iter, temp = y2,
        tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

    cur = obj.fun(res1$par)

    print(test.obj.fun.new(res1$par))
    
    if(cur == 100) return(res1$par[-length(res1$par)])

    while ((cur - old) > 1e-06) {
        res1 <- optim(res1$par, obj.fun, swap.stage2, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        res1 <- optim(res1$par, obj.fun, swap.stage1, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        res1 <- optim(res1$par, obj.fun, swap.stage3, method = "SANN", control = list(maxit = iter, temp = y2,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        old = cur
        cur = obj.fun(res1$par)

        print(test.obj.fun.new(res1$par))


        if(cur == 100) return(res1$par[-length(res1$par)])
    }

    c(nTrt, bRep, tRep, nPlot)

    inter.stage1 = exp(log(y2.stage1/y2low.stage1)/8)
    inter.stage2 = exp(log(y2.stage2/y2low.stage2)/8)
    inter = exp(log(y2/y2low)/8)

    i = 1
    unstable = FALSE

    # for (i in 1:8) {
    while (i < 9 || unstable) {
        unstable = FALSE
        y2.stage2 = y2.stage2/inter.stage2
        y2.stage1 = y2.stage1/inter.stage1
        y2 = y2/inter
        
        old = obj.fun(res1$par)

        cat("Level: ", i + 2, ", Current temp: ", y2.stage2, ",", y2.stage1, ",", y2, "\n")
        res1 <- optim(res1$par, obj.fun, swap.stage2, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        res1 <- optim(res1$par, obj.fun, swap.stage1, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
        
        res1 <- optim(res1$par, obj.fun, swap.stage3, method = "SANN", control = list(maxit = iter, temp = y2,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        cur = obj.fun(res1$par)
          print(test.obj.fun.new(res1$par))

        if(cur == 100) return(res1$par[-length(res1$par)])

        while ((cur - old) > 1e-06) {
            unstable = TRUE
            res1 <- optim(res1$par, obj.fun, swap.stage2, method = "SANN", control = list(maxit = iter,
                temp = y2.stage2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

            res1 <- optim(res1$par, obj.fun, swap.stage1, method = "SANN", control = list(maxit = iter,
                temp = y2.stage1, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

            res1 <- optim(res1$par, obj.fun, swap.stage3, method = "SANN", control = list(maxit = iter, temp = y2,
              tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

            old = cur
            cur = obj.fun(res1$par)
            print(test.obj.fun.new(res1$par))

            if(cur == 100) return(res1$par[-length(res1$par)])

        }


        i = i + 1
    }

    print(c(nTrt, bRep, tRep, nPlot))

    return(res1$par[-length(res1$par)])
}

set.seed(1001)
bestInd = threeStageCRD(init =  sample(1:n), iter = 3000, obj.fun = obj.fun.CRD, 
test.obj.fun.new= test.obj.fun.CRD)

############################################################################################

test.obj.fun.CRD = function(ind) {

        
    #inforamtion in within runs and tags stratum while the cage information is eliminated
   
       Rep = Z1.rep
       newZ1 = Z1.mat[ind[-length(ind)], ]
       PP = projMat(newZ1) %*% blk.proj  
 
        Rep = trt.rep
       newX = X.trt[ind[-length(ind)], ]
       info.mat = t(newX) %*% blk.proj %*% newX
       e.va = eigen(info.mat, only.values = TRUE)$va
       can.eff = e.va[-which(e.va < 1e-07)]/Rep
       len = length(can.eff)

    
    return(list(aveEff = 100/mean(1/can.eff), res = tr(PP) - len))
}

################################################################################
#Initialise the temperature

initialTemp = function(iter, newInit, swap, obj.fun){
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  U = numeric(iter)
  temp = matrix(0, nrow = iter, ncol = length(newInit))
  temp[1,] = newInit
  U[1] =   obj.fun(temp[1,])
  tol = 10
  check.U = numeric(tol) 
  
  for(i in 2:iter){
    setTxtProgressBar(pb, i)
    temp[i,] = swap(temp[i-1,])
    U[i] = obj.fun(temp[i,])
    
    check.U = c(check.U, U[i])
    check.U = check.U[-1]
    
    if(all(outer(check.U,check.U,"=="))) temp[i,] = newInit
  }
  close(pb)
  
  #if(any(U == 0)) U = U[-which(U ==0)]
  
  #y2 = max(U) - min(U) 
  
  #y2low = min(range((abs(diff(U))[which(abs(diff(U))>1e-5)])))
  if(all(U==0) || all(diff(U)<1e-8)){
    return(list(y2 = 100, y2low = 1e-5, newInit = newInit))
  
  }else{
    xxx = max(U, na.rm = TRUE) - range(U[which(U < (max(U, na.rm = TRUE)-1e-10))], na.rm = TRUE)
    
    if(isTRUE(all.equal(xxx[1], xxx[2]))) xxx[2] = 1e-5
     
    newInit = temp[which(U == max(U))[1],]
    return(list(y2 = max(xxx), y2low = min(xxx), newInit =newInit))
  }
}   


#######################################################################################################
#Final function

optThreeStage =  function(init, iter = 10000, obj.fun, test.obj.fun) {
    print(c(nTrt, bRep, tRep, nPlot))
   
    newInit = c(init, init[1])

    old = obj.fun(newInit)

    cat("Level: 1, Finding the temperatures for both stages\n")

    temp = initialTemp(iter, newInit, swap.stage2.new, obj.fun)

    res.stage2 = temp$newInit
    y2.stage2 = temp$y2
    y2low.stage2 = temp$y2low


    temp = initialTemp(iter, newInit, swap.stage1.new, obj.fun)

    res.stage1 = temp$newInit
    y2.stage1 = temp$y2
    y2low.stage1 = temp$y2low


    temp = initialTemp(iter, newInit,swap.stage3.new, obj.fun)

    res = temp$newInit
    y2 = temp$y2
    y2low = temp$y2low
    
    cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
    cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
    cat("Stage 0 temperature range: ", y2, " to ", y2low, "\n")


    test.newInit = apply(rbind(res,res.stage1,res.stage2),1, obj.fun)
    newInit = rbind(res,res.stage1,res.stage2)[which(test.newInit == max(test.newInit))[1],]
    
    cur = obj.fun(newInit)
    print(test.obj.fun(newInit))
    
    while ((cur - old) > 1e-06) {
        
        temp = initialTemp(iter, newInit, swap.stage2.new, obj.fun)
    
        res.stage2 = temp$newInit
        y2.stage2 = temp$y2
        y2low.stage2 = temp$y2low
    
    
        temp = initialTemp(iter, newInit, swap.stage1.new, obj.fun)
    
        res.stage1 = temp$newInit
        y2.stage1 = temp$y2
        y2low.stage1 = temp$y2low
    
    
        temp = initialTemp(iter, newInit, swap.stage3.new, obj.fun)
    
        res = temp$newInit
        y2 = temp$y2
        y2low = temp$y2low
           
        test.newInit = apply(rbind(res,res.stage1,res.stage2),1, obj.fun)
        newInit = rbind(res,res.stage1,res.stage2)[which(test.newInit == max(test.newInit))[1],]

        old = cur
        
        cur = obj.fun(newInit)
        print(test.obj.fun(newInit))     
   
        cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
        cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
        cat("Stage 0 temperature range: ", y2, " to ", y2low, "\n")
    
     }
      

    old = obj.fun(newInit)
  
    cat("Level: 2, Current temp: ", y2.stage2, ",", y2.stage1, ",", y2, "\n")

    res <- optim(newInit, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
        tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

    res1 <- optim(res$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
        tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

    res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, temp = y2,
        tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

    cur = obj.fun(res1$par)

    print(test.obj.fun(res1$par))
    
    if(cur == 100) return(res1$par[-length(res1$par)])

    while ((cur - old) > 1e-06) {
        res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        res1 <- optim(res1$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, temp = y2,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        old = cur
        cur = obj.fun(res1$par)

        print(test.obj.fun(res1$par))


        if(cur == 100) return(res1$par[-length(res1$par)])
    }

    c(nTrt, bRep, tRep, nPlot)

    inter.stage1 = exp(log(y2.stage1/y2low.stage1)/8)
    inter.stage2 = exp(log(y2.stage2/y2low.stage2)/8)
    inter = exp(log(y2/y2low)/8)

    i = 1
    unstable = FALSE

    # for (i in 1:8) {
    while (i < 9 || unstable) {
        unstable = FALSE
        y2.stage2 = y2.stage2/inter.stage2
        y2.stage1 = y2.stage1/inter.stage1
        y2 = y2/inter
        
        old = obj.fun(res1$par)

        cat("Level: ", i + 2, ", Current temp: ", y2.stage2, ",", y2.stage1, ",", y2, "\n")
        res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        res1 <- optim(res1$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))
        
        res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, temp = y2,
            tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

        cur = obj.fun(res1$par)
          print(test.obj.fun(res1$par))

        if(cur == 100) return(res1$par[-length(res1$par)])

        while ((cur - old) > 1e-06) {
            unstable = TRUE
            res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter,
                temp = y2.stage2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

            res1 <- optim(res1$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter,
                temp = y2.stage1, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

            res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, temp = y2,
              tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

            old = cur
            cur = obj.fun(res1$par)
            print(test.obj.fun(res1$par))

            if(cur == 100) return(res1$par[-length(res1$par)])

        }


        i = i + 1
    }

    print(c(nTrt, bRep, tRep, nPlot))

    return(res1$par[-length(res1$par)])
}

                               
################################################################################
#Compare Two-stages swapping method versus standard swapping method. 
init = 1:n #sample(1:n)

obj.fun.new(c(1:n,1))

old.time = proc.time()
#bestInd = newInit[-length(newInit)]
#bestInd = res$par[-length(res$par)]
set.seed(527)
bestInd = optThreeStage(init =  1:n, iter = 10000, obj.fun = obj.fun.CRD, test.obj.fun = test.obj.fun.CRD)
proc.time() - old.time

#bestInd = 1:n

test(X.trt = Z1.mat[bestInd,], blk.proj, Rep=Z1.rep)$ave.eff
test(X.trt = Z1.mat[bestInd,], blk.proj, Rep=Z1.rep)$can.eff

test(X.trt = Z1.mat[bestInd,], mI(n) - Pb, Rep=Z1.rep)$ave.eff
test(X.trt = Z1.mat[bestInd,], mI(n) - Pb, Rep=Z1.rep)$can.eff

test(X.trt = X.trt[bestInd,], mI(n) - Pb, Rep=trt.rep)$can.eff
test(X.trt = X.trt[bestInd,], blk.proj, Rep=trt.rep)$can.eff

test(X.trt = X.trt[bestInd,], Pb - mK(n), Rep=trt.rep)$can.eff
test(X.trt = X.trt[bestInd,], Pb1 - mK(n), Rep=trt.rep)$can.eff


new.Z1.mat=  Z1.mat[bestInd,]
#new.Z1.mat=  Z1.mat

################################################################################
#Set the design from the search 

colnames(new.Z1.mat) = sort(levels(interaction(multiLetters(1:nZ1))))
new.Z1.des = apply(new.Z1.mat, 1,  function(x)
colnames(new.Z1.mat)[which(as.logical(x))])
new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))

design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

design.df = cbind(design.df,
            Ani = t(new.Z1.des),
            Pro = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                          as.character(phase1DesignEX1$Ani)),]$Pro,
            Sam = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                          as.character(phase1DesignEX1$Ani)),]$Sam)

################################################################################
#Theortical ANOVA table 

summaryAovOnePhase(design.df, blk.str = "Ani", 
trt.str = "Pro")


summaryAovOnePhase(design.df, blk.str = "Ani", 
trt.str = "Trt")

summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", 
trt.str = "Tag * Pro")

summaryAovOnePhase(design.df, blk.str = "Ani", 
trt.str = "Trt", latex = TRUE, fixed.names= "\\tau")

summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", 
trt.str = "Tag + Trt", latex = TRUE)


summaryAovOnePhase(design.df, blk.str = "Run", 
trt.str = "Tag")

summaryAovTwoPhase(design.df,  blk.str2 = "Tag", blk.str1 = "Ani",
trt.str = "Trt")

summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Trt")

summaryAovTwoPhase(design.df,  blk.str2 = "Run + Tag", blk.str1 = "Ani",
trt.str = "Trt")

summaryAovOnePhase(design.df,  blk.str = "Run + Tag", trt.str = "Trt")

summaryAovOnePhase(design.df,  blk.str = "Run + Tag", trt.str = "Ani")

matrix(design.df$Ani, nrow = nBlk, ncol = nPlot, byrow = TRUE)
(N = with(design.df, table(Ani, Run)))
N %*% t(N)

(N = with(design.df, table(Ani, Tag)))
N %*% t(N)

matrix(tolower(design.df$Trt), nrow = nBlk, ncol = nPlot, byrow = TRUE)
(N = with(design.df, table(Trt, Run)))
N %*% t(N)

(N = with(design.df, table(Trt, Tag)))
N %*% t(N)



fileName = paste("design_", nTrt, "trt_", nAni, "ani_", nPlot, "tag_", n/nPlot, 
"run.Rdata", sep = "")

fileName

save(design.df, file = fileName)


