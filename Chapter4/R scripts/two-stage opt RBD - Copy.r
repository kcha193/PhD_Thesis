nTrt = 8 
bRep = 10
nCag = 10
tRep = 2
n = nTrt * bRep * tRep
nPlot = 4
nBlk = n/nPlot
trt.rep = n/nTrt

nAni = nTrt * bRep
nBlk = n/nPlot
nZ1 =  nAni
Z1.rep = n/nZ1

c(nTrt, nAni, nPlot, n/nPlot)
  n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt

    nAni = nTrt * bRep

    if (!is.wholenumber(n/nPlot)) {
        stop("The samples from Phase 1 experiment can not be fitted to the Phase 2 experiment!")
    }


    cat("Design parameters:\n")
    cat("Trt:", nTrt, "Ani:", nAni, "Tag:", nPlot, "Run:", nBlk, "\n")

    phase1DesignEX1 <- local({

        Cag = rep(1:nCag, each = (nTrt * bRep)/nCag)
        Ani = 1:nAni
        Trt = 1:nTrt

        data.frame(cbind(Cag, Ani, Trt))
    })

    phase1DesignEX1$Ani = as.factor(newLETTERS[phase1DesignEX1$Ani])
    phase1DesignEX1$Trt = as.factor(multiLetters(phase1DesignEX1$Trt))
    phase1DesignEX1$Cag = as.factor(letters[phase1DesignEX1$Cag])

    phase1DesignEX1

#########################################
#orthogonal projectors

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
#########################################




Z1.des = initialRBD(nTrt = nTrt, bRep = bRep, nCag = nCag, tRep = tRep, nPlot = nPlot)


############################################################

Z1.mat = matrix(0, ncol = nAni, nrow = nAni * Z1.rep)
Z1.mat[cbind(1:(nAni * Z1.rep), Z1.des)] = 1
Z1.rep = nrow(Z1.mat)/ncol(Z1.mat)



# Parameter's of treatment
trt.rep = n/nTrt

trt.des = as.numeric(phase1DesignEX1$Trt[match(newLETTERS[Z1.des], phase1DesignEX1$Ani)])
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
    #Two checks for omitting the swap of two indtical animals and animals in the different block
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

obj.fun.old = function(ind) {
     
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

resDF = tr(projMat(Z1.mat)) - (nCag - 1) - (nTrt -1) - 1
C.cage  = (mI(nCag) - mK(nCag)) %x%  mK(nAni/nCag) 
C.ani =   mI(nCag)%x%  (mI(nAni/nCag) - mK(nAni/nCag))
cage.Rep = n/nCag

C.ani =   mI(nCag)%x%  (mI(nAni/nCag) - mK(nAni/nCag))
ani.Rep = n/nAni
trt.Rep = trt.rep
  
#Objective functions
obj.fun.old1 = function(ind) {

 info.mat = matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)
    
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    can.eff = e.va[-which(e.va < 1e-07)]/cage.Rep 
    len.cage = length(can.eff)
    aveEff2 = 100/mean(1/can.eff)
    if(is.na(aveEff2)) aveEff1 = 0
 
 
    if(any(abs(info.mat) > 1e-7)){     #check for any information in the cage stratum
      PP = blk.proj - matMulti1(blk.proj, ginv( matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)), C.cage,  Z1.mat[ind[-length(ind)],]) 
    } else{       #if there is no information in the cage stratum then no decomposition is required. If decomposition is performed here, 
                  # some very werid thing is going to happen. This is because taking the inverse of a very small number will become huge. 
                  # Hence, large information of the within block within runs can be extracted. 
      PP = blk.proj
    }
    
    #Maximise the amount of animal information in within runs and tags
    info.mat = matMulti(PP , Z1.mat[ind[-length(ind)],], C.ani)

    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    can.eff = e.va[-which(e.va < 1e-07)]/ani.Rep 

    aveEff1 = 1/mean(1/can.eff)
    if(is.na(aveEff1)) aveEff1 = 0
 
      
   #Maximise the amount of treatment information in within runs and tags

    PP =  matMulti1(blk.proj, ginv(info.mat), C.ani, Z1.mat[ind[-length(ind)],]) 
           
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*%  PP %*% newX
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    can.eff = e.va[-which(e.va < 1e-07)]/trt.Rep
    len = length(can.eff)
        
    aveEff = (3/4) * aveEff1 + 1/4 * ((len + 1/mean(1/can.eff))/nTrt)        

    return(aveEff)
}


###################################################################################
resDF = tr(projMat(Z1.mat)) - (nCag - 1) - (nTrt -1) - 1
C.cage  = (mI(nCag) - mK(nCag)) %x%  mK(nAni/nCag) 
C.ani =   mI(nCag)%x%  (mI(nAni/nCag) - mK(nAni/nCag))
cage.Rep = n/nCag

C.ani =   mI(nCag)%x%  (mI(nAni/nCag) - mK(nAni/nCag))
ani.Rep = n/nAni
trt.Rep = trt.rep

obj.fun.new1 = function(ind) {
        
    #inforamtion in within runs and tags stratum while the cage information is eliminated
    info.mat = matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage) 
 
    if(any(abs(info.mat) > 1e-7)){     #check for any information in the cage stratum
      PP = blk.proj - matMulti1(blk.proj, ginv( matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)), C.cage,  Z1.mat[ind[-length(ind)],]) 
    } else{       #if there is no information in the cage stratum then no decomposition is required. If decomposition is performed here, 
                  # some very werid thing is going to happen. This is because taking the inverse of a very small number will become huge. 
                  # Hence, large information of the within block within runs can be extracted. 
      PP = blk.proj
    }
    
    #Maximise the amount of animal information in within runs and tags
    info.mat = matMulti(PP , Z1.mat[ind[-length(ind)],], C.ani)
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    can.eff = e.va[-which(e.va < 1e-07)]/ani.Rep 

    aveEff1 = 1/mean(1/can.eff)
    
    if(is.na(aveEff1)) aveEff1 = 0
 
      
   #Maximise the amount of treatment information in within runs and tags

    PP =  matMulti1(blk.proj, ginv(info.mat), C.ani, Z1.mat[ind[-length(ind)],]) 
           
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*%  PP %*% newX
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    can.eff = e.va[-which(e.va < 1e-07)]/trt.Rep
    len = length(can.eff)
               
    aveEff = (((tr(PP) - len)/resDF) * (1/resDF) +#residual DF 
             aveEff1 * (resDF - 1)/resDF) * 3/4 +             #Animal information
             1/4 * ((len + 1/mean(1/can.eff))/nTrt)  #Treatment information

    return(aveEff)
}

###################################################################################


   #aveEff =  ((tr(PP) - len)/resDF) *   (100/9) +#residual DF
    #         aveEff1 * (5/9) +                               #Animal information
    #         ((1/mean(1/can.eff)) *(1/nTrt) + len/(nTrt-1) * ((nTrt-1)/nTrt)) * 100/3 #Treatment information
    
    #aveEff =  ((tr(PP) - len)/resDF) *   (100/(3*nTrt)) +#residual DF
    #         aveEff1 * (2/3 - 1/(3*nTrt)) +                               #Animal information
    #         ((1/mean(1/can.eff)) *(1/nTrt) + len/(nTrt-1) * ((nTrt-1)/nTrt)) * 100/3 #Treatment information


    #aveEff =  (((tr(PP) - len)/resDF) * (100/resDF) + #residual DF
    #         aveEff1 * (resDF - 1)/resDF) * 2/3 +                               #Animal information
    #         ((1/mean(1/can.eff)) *(1/nTrt) +
    #         len/(nTrt-1) * ((nTrt-1)/nTrt)) * 100/3 #Treatment information
 
########################################################################################################### 
 
resDF = tr(projMat(Z1.mat)) - (nCag - 1) - (nTrt -1) - 1
C.cage  = (mI(nCag) - mK(nCag)) %x%  mK(nAni/nCag) 
C.ani =   mI(nCag)%x%  (mI(nAni/nCag) - mK(nAni/nCag))
cage.Rep = n/nCag

weight = resDF #ceiling(sqrt(resDF))

C.ani =   mI(nCag)%x%  (mI(nAni/nCag) - mK(nAni/nCag))
ani.Rep = n/nAni
trt.Rep = trt.rep

try.obj.fun.new = function(ind) {

    #inforamtion in within runs and tags stratum while the cage information is eliminated
    info.mat = matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)
    PP = matMulti1(blk.proj, ginv(info.mat), C.cage,  Z1.mat[ind[-length(ind)],])
      
    if(any(abs(info.mat) > 1e-7)){     
    #check for any information in the cage stratum
      PP = blk.proj - PP 
    } else{      
     #if there is no information in the cage stratum then no decomposition is required. 
     #If decomposition is performed here, 
     # some very werid thing is going to happen. 
     #This is because taking the inverse of a very small number will become huge. 
     # Hence, large information of the within block within runs can be extracted. 
      PP = blk.proj
    }
    
    #Maximise the amount of animal information in within runs and tags
    info.mat = matMulti(blk.proj , Z1.mat[ind[-length(ind)],], C.ani)
   
    #Maximise the amount of animal information in within runs and tags
  
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    can.eff = e.va[-which(e.va < 1e-07)]/ani.Rep 

    aveEff1 = 1/mean(1/can.eff)
     
    if(is.na(aveEff1)) aveEff1 = 0
   
    if(!isTRUE(all.equal(aveEff1, 1))) return(0)
   
    #Maximise the amount of treatment information in within runs and tags
    PP =  matMulti1(blk.proj, ginv(info.mat), C.ani, Z1.mat[ind[-length(ind)],]) 

    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*%  PP %*% newX
    #len = length(can.eff)

    #len = qr(info.mat)$rank 
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
 
    can.eff = e.va[-which(e.va < 1e-07)]/trt.Rep
    len = length(can.eff)
  
    if(!isTRUE(all.equal(len/(nTrt-1), 1))) return(0)
 
    aveEff =  tr(PP) - len 
    
    return(aveEff)
}


obj.fun.new = function(ind) {

  ans = try(try.obj.fun.new(ind), silent = TRUE)
  
  ifelse(class(ans) == "try-error", 0, ans)
}

test.obj.fun.new = function(ind) {

        
    #inforamtion in within runs and tags stratum while the cage information is eliminated
    info.mat = matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)
    
    if(any(abs(info.mat) > 1e-7)){     #check for any information in the cage stratum
      PP = blk.proj %*%  blk.proj - matMulti1(blk.proj, ginv( matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)), C.cage,  Z1.mat[ind[-length(ind)],]) 
    } else{       #if there is no information in the cage stratum then no decomposition is required. If decomposition is performed here, 
                  # some very werid thing is going to happen. This is because taking the inverse of a very small number will become huge. 
                  # Hence, large information of the within block within runs can be extracted. 
      PP = blk.proj
    }
    
    #Maximise the amount of animal information in within runs and tags
    info.mat = matMulti(PP , Z1.mat[ind[-length(ind)],], C.ani)
      
      
   #Maximise the amount of treatment information in within runs and tags

    PP =  matMulti1(blk.proj, ginv(info.mat), C.ani, Z1.mat[ind[-length(ind)],]) 
           
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*%  PP %*% newX
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    #e.va = Re(svd(info.mat)$d)
    can.eff = e.va[-which(e.va < 1e-07)]/trt.Rep
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
  
  if(all(U==0) || all(diff(U)<1e-8)){
    return(list(y2 = 100, y2low = 1e-5, newInit = newInit))
  
  }else{
    xxx = max(U, na.rm = TRUE) - range(U[which(U < (max(U, na.rm = TRUE)-1e-10))], na.rm = TRUE)
    
    if(isTRUE(all.equal(xxx[1], xxx[2]))) xxx[2] = 1e-5
     
    newInit = temp[which(U == max(U, na.rm = TRUE))[1],]
    return(list(y2 = max(xxx), y2low = min(xxx), newInit =newInit))
  }
}   

 
      
################################################################################
#Two-stages swapping method

threeStageRBD = function(init, iter = 50000, obj.fun = obj.fun.new) {
    print(c(nTrt, bRep, nCag, tRep, nPlot))
   
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


    temp = initialTemp(iter, newInit, swap.stage3.new, obj.fun)

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
        print(test.obj.fun.new(newInit))     
   
        cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
        cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
        cat("Stage 0 temperature range: ", y2, " to ", y2low, "\n")
    
     }
      

    old = obj.fun(newInit)
  
    cat("Level: 2, Current temp: ", y2.stage2, ",", y2.stage1, ",", y2, "\n")

    res <- optim(newInit, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
        tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

    res1 <- optim(res$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
        tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

    res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, temp = y2,
        tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

    cur = obj.fun(res1$par)

    print(test.obj.fun.new(res1$par))
    
    if(cur == 1) return(res1$par[-length(res1$par)])

    while ((cur - old) > 1e-06) {
        res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

        res1 <- optim(res1$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

        res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, temp = y2,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

        old = cur
        cur = obj.fun(res1$par)

        print(test.obj.fun.new(res1$par))


        if(cur == 1) return(res1$par[-length(res1$par)])
    }


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
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

        res1 <- optim(res1$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))
        
        res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, temp = y2,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

        cur = obj.fun(res1$par)
          print(test.obj.fun.new(res1$par))

        if(cur == 1) return(res1$par[-length(res1$par)])

        while ((cur - old) > 1e-06) {
            unstable = TRUE
            res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter,
                temp = y2.stage2, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

            res1 <- optim(res1$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter,
                temp = y2.stage1, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

            res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, temp = y2,
              tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

            old = cur
            cur = obj.fun(res1$par)
            print(test.obj.fun.new(res1$par))

            if(cur == 1) return(res1$par[-length(res1$par)])

        }


        i = i + 1
    }

    print(c(nTrt, bRep, nCag, tRep, nPlot))

    return(res1$par[-length(res1$par)])
}

#res1 <- optim(res1$par, obj.fun.new, swap.new, method = "SANN", control = list(maxit = iter,
#        temp = y2, tmax = 1000, trace = TRUE, REPORT = 10, fnscale = -1))

c(nTrt, bRep, nCag, tRep, nPlot)
#bestInd = res1$par[-length(res1$par)]

################################################################################
#Objective function section comparesion 1
#With Example nTrt = 3 bRep = 4 nCag = 4 tRep = 2

print(obj.fun.old(c(1:n,1)))

bestInd1 = threeStageRBD(init =  1:n, iter = 1000, obj.fun = obj.fun.old)

quickTest(bestInd1)

print(obj.fun.old1(c(1:n,1)))

bestInd2 = threeStageRBD(init =  1:n, iter = 1000, obj.fun = obj.fun.old1)

quickTest(bestInd2)
################################################################################
#Objective function section comparesion 2
#With Example nTrt = 4 bRep = 2 nCag = 2 tRep = 2

print(obj.fun.old1(c(1:n,1)))

bestInd1 = threeStageRBD(init =  1:n, iter = 1000, obj.fun = obj.fun.old1)

quickTest(bestInd1)

print(obj.fun.new1(c(1:n,1)))

bestInd2 = threeStageRBD(init =  1:n, iter = 1000, obj.fun = obj.fun.new1)

quickTest(bestInd2)
################################################################################
#Objective function section comparesion 3
#With Example nTrt = 4 bRep = 2 nCag = 2 tRep = 2  nPlot = 4
print(obj.fun.new1(c(1:n,1)))

bestInd1 = threeStageRBD(init =  1:n, iter = 1000, obj.fun = obj.fun.new1)

quickTest(bestInd1)

print(obj.fun.new(c(1:n,1)))

bestInd1 = threeStageRBD(init =  1:n, iter = 1000, obj.fun = obj.fun.new)


print(test.obj.fun.new(c(bestInd1,1)))

if(!isTRUE(all.equal(obj.fun.new(c(1:n,1)),0)) &&
  isTRUE(all.equal(test.obj.fun.new(c(1:n,1))$res, test.obj.fun.new(c(bestInd1,1))$res)) && 
  (test.obj.fun.new(c(1:n,1))$aveEff > test.obj.fun.new(c(bestInd1,1))$aveEff)){
  bestInd1 = 1:n
}

newResDF = test.obj.fun.new(c(bestInd1,1))$res
# newResDF =  newResDF -1
################################################################################

try.obj.fun.new1 = function(ind) {

    #inforamtion in within runs and tags stratum while the cage information is eliminated
    info.mat = matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)
    PP = matMulti1(blk.proj, ginv(info.mat), C.cage,  Z1.mat[ind[-length(ind)],])
         
    if(any(abs(info.mat) > 1e-7)){     #check for any information in the cage stratum
      PP = blk.proj - PP 
    } else{     
      PP = blk.proj
    }
    
    #Maximise the amount of animal information in within runs and tags
    info.mat = matMulti(PP , Z1.mat[ind[-length(ind)],], C.ani)
  
    PP =  matMulti1(blk.proj, ginv(info.mat), C.ani, Z1.mat[ind[-length(ind)],]) 
     
    aveEff1 = tr(info.mat/ani.Rep)/tr(PP)
     
    if(is.na(aveEff1)) aveEff1 = 0
   
    if(!isTRUE(all.equal(aveEff1, 1))) return(0)
            
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*%  PP %*% newX
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    #e.va = Re(svd(info.mat)$d)
    can.eff = e.va[-which(e.va < 1e-07)]/trt.Rep
    len = length(can.eff)

    if(!isTRUE(all.equal(len/(nTrt-1), 1))) return(0)
    if(!isTRUE(all.equal((tr(PP) - len), newResDF))) return(0)  
   
    aveEff = 100/(mean(1/can.eff))
    #aveEff = tr(info.mat %*% info.mat)
    
    return(aveEff)
}


obj.fun.new1 = function(ind) {

  ans = try(try.obj.fun.new1(ind), silent = TRUE)
  
  ifelse(class(ans) == "try-error", 0, ans)
}

obj.fun.new1(c(bestInd1,1))
 
old = proc.time()
#bestInd1 = 1:n 
bestInd2 = threeStageRBD(init = bestInd1, iter = 1000, obj.fun = obj.fun.new1)
proc.time() - old 

quickTest(bestInd2)

##########################################################################################
#SA comparison with the initial design 
#With Example nTrt = 7 bRep = 4 nCag = 2 tRep = 2  nPlot = 4


print(test.obj.fun.new(c(1:n,1)))

obj.fun.new(c(1:n,1))

bestInd1 = threeStageRBD(init =  1:n, iter = 1000, obj.fun = obj.fun.new)

print(test.obj.fun.new(c(bestInd1,1)))

if(!isTRUE(all.equal(obj.fun.new(c(1:n,1)),0)) &&
  isTRUE(all.equal(test.obj.fun.new(c(1:n,1))$res, test.obj.fun.new(c(bestInd1,1))$res)) && 
  (test.obj.fun.new(c(1:n,1))$aveEff > test.obj.fun.new(c(bestInd1,1))$aveEff)){
  bestInd1 = 1:n
}

newResDF = test.obj.fun.new(c(bestInd1,1))$res
# newResDF =  newResDF -1
################################################################################

try.obj.fun.new1 = function(ind) {

    #inforamtion in within runs and tags stratum while the cage information is eliminated
    info.mat = matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)
    PP = matMulti1(blk.proj, ginv(info.mat), C.cage,  Z1.mat[ind[-length(ind)],])
         
    if(any(abs(info.mat) > 1e-7)){     #check for any information in the cage stratum
      PP = blk.proj - PP 
    } else{     
      PP = blk.proj
    }
    
    #Maximise the amount of animal information in within runs and tags
    info.mat = matMulti(PP , Z1.mat[ind[-length(ind)],], C.ani)
  
    PP =  matMulti1(blk.proj, ginv(info.mat), C.ani, Z1.mat[ind[-length(ind)],]) 
     
    aveEff1 = tr(info.mat/ani.Rep)/tr(PP)
     
    if(is.na(aveEff1)) aveEff1 = 0
   
    if(!isTRUE(all.equal(aveEff1, 1))) return(0)
            
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*%  PP %*% newX
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    #e.va = Re(svd(info.mat)$d)
    can.eff = e.va[-which(e.va < 1e-07)]/trt.Rep
    len = length(can.eff)

    if(!isTRUE(all.equal(len/(nTrt-1), 1))) return(0)
    if(!isTRUE(all.equal((tr(PP) - len), newResDF))) return(0)  
   
    aveEff = 100/(mean(1/can.eff))
    #aveEff = tr(info.mat %*% info.mat)
    
    return(aveEff)
}


obj.fun.new1 = function(ind) {

  ans = try(try.obj.fun.new1(ind), silent = TRUE)
  
  ifelse(class(ans) == "try-error", 0, ans)
}

obj.fun.new1(c(bestInd1,1))
 
old = proc.time()
#bestInd1 = 1:n 
bestInd2 = twoStageRBD(init = bestInd1, iter = 10000, obj.fun = obj.fun.new1)
proc.time() - old 

 old = proc.time()

bestInd2 = twoStageRBD(init = bestInd2, iter = 10000, obj.fun = obj.fun.new1)
 proc.time() - old 

print(obj.fun.new1(c(bestInd2,1)))

bestInd = bestInd2
#################################################################################
info.mat =  matMulti(( Pb - mK(n)) , Z1.mat[bestInd,], C.cage)

info.mat =  matMulti((mI(n) - Pb) , Z1.mat[bestInd,], C.cage)

if(any(abs(info.mat) > 1e-7)){
  PP = (mI(n) - Pb) - matMulti1((mI(n) - Pb), ginv( matMulti((mI(n) - Pb), 
      Z1.mat[bestInd,], C.cage)), C.cage, Z1.mat[bestInd,]) 
  
  PP1 = matMulti1(PP, ginv( matMulti(PP , Z1.mat[bestInd,], C.ani)), 
      C.ani, Z1.mat[bestInd,]) 
      
} else{
  PP1 = matMulti1((mI(n) - Pb), ginv( matMulti((mI(n) - Pb), 
      Z1.mat[bestInd,], C.ani)), C.ani, Z1.mat[bestInd,]) 
} 
 
test(X.trt = mI(n), blk.proj = PP1, Rep = 1)[c("nCan", "can.eff")]


info.mat =  matMulti(blk.proj , Z1.mat[bestInd,], C.cage)

if(any(abs(info.mat) > 1e-7)){
PP = blk.proj - matMulti1(blk.proj, ginv( matMulti(blk.proj , Z1.mat[bestInd,], C.cage)), C.cage, Z1.mat[bestInd,]) 

PP1 = matMulti1(PP, ginv( matMulti(PP , Z1.mat[bestInd,], C.ani)), C.ani, Z1.mat[bestInd,]) 
} else{
PP1 = matMulti1( blk.proj, ginv( matMulti( blk.proj , Z1.mat[bestInd,], C.ani)), C.ani, Z1.mat[bestInd,]) 

} 

test(X.trt = X.trt[bestInd,], PP1, Rep=trt.rep)


fractions(test(X.trt = X.trt[bestInd,], PP1, Rep=trt.rep)$can)

new.Z1.mat=  Z1.mat[bestInd,]
#new.Z1.mat=  Z1.mat

################################################################################
#Set the design from the search

colnames(new.Z1.mat) = sort(levels(interaction(newLETTERS[1:nZ1])))
new.Z1.des = apply(new.Z1.mat, 1,  function(x)
colnames(new.Z1.mat)[which(as.logical(x))])
new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))

design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

design.df = cbind(design.df,
            Cag = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                          as.character(phase1DesignEX1$Ani)),]$Cag,
            Ani = t(new.Z1.des),
            Trt = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                          as.character(phase1DesignEX1$Ani)),]$Trt)

summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt")

levels(design.df$Ani) = rep(1:(nlevels(design.df$Ani)/nCag), nCag)

summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")

################################################################################
c(nTrt, bRep, nCag, tRep, nPlot)

fileName = paste("RBDdesign_", nTrt, "trt_", nAni, "ani_", nCag,"blk_",  nPlot, "tag_", 
n/nPlot, "run.Rdata", sep = "")

fileName

save(design.df, file = fileName)

################################################################################


summary.aov.onePhase(design.df,  blk.str = "Run", trt.str = "Cag/Ani")
summary.aov.onePhase(design.df,  blk.str = "Run + Tag", trt.str = "Cag")

summary.aov.onePhase(design.df,  blk.str = "Run", trt.str = "Ani")

summary.aov.onePhase(design.df,  blk.str = "Cag/Ani", trt.str = "Trt")

summary.aov.onePhase(design.df,  blk.str = "Cag", trt.str = "Ani")

matrix(design.df$Cag, nrow = nBlk, ncol = nPlot, byrow = TRUE)

matrix(design.df$Ani, nrow = nBlk, ncol = nPlot, byrow = TRUE)

matrix(design.df$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE)

matrix(paste(toupper(design.df$Cag), design.df$Ani, sep = ""), nrow = nBlk, ncol = nPlot, byrow = TRUE)


matrix(design.df$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE)

N = with(design.df, table(Trt,Run))
N %*% t(N)

N = with(design.df, table(Trt,Cag))
N %*% t(N)

(N = with(design.df, table(Cag,Run)))
N %*% t(N)

(N = with(design.df, table(Ani,Run)))
N %*% t(N)


N =with(design.df, table(paste(Cag, Ani, sep = ""), Run))
CM = N %*% t(N)

CM %*% CM %*% CM %*% CM%*% CM %*% CM %*% CM

with(design.df, table(Trt, Cag, Run))

N = with(design.df, table(Ani,Run))
N %*% t(N)

 summary.aov.twoPhase(design.df,  blk.str2 = "Run+ Tag", blk.str1 = "Cag/Ani",
trt.str = "Trt")

summary.aov.twoPhase(design.df,  blk.str2 = "Run+ Tag", blk.str1 = "Ani",
trt.str = "Trt")

summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = " Trt")

summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag",
trt.str = "Tag + Trt")


summary.aov.twoPhase(design.df,  blk.str2 = "Run + Tag", blk.str1 = "Cag/Ani",
trt.str = "Trt")



des1 = design.df

                                
des1$Cag = factor(rep(letters[1:3],each = 12))
 des1$Ani = factor(c( 1,2,3,4,
                      2,1,4,3,
                      5,5,6,6,
                      3,6,1,5,
                      6,3,5,1,
                      4,4,2,2,
                      5,2,6,4,
                      2,5,4,6,
                      1,1,3,3))
des1$Trt = factor(letters[
            c(1,2,3,4,
              2,1,4,3,
              5,5,6,6,
              3,6,1,5,
              6,3,5,1,
              4,4,2,2,
              5,2,6,4,
              2,5,4,6,
              1,1,3,3)])

summary.aov.twoPhase(des1,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt")


N = with(des1, table(Trt,Run))
N %*% t(N)


matrix(des1$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE)



trt.contr = test(X.trt = X.trt[bestInd,],  PP1, Rep=trt.rep)$e.vec

trt.contr1 = trt.contr[,1][as.numeric(design.df$Trt)]
trt.contr2 = trt.contr[,2][as.numeric(design.df$Trt)]
trt.contr3 = trt.contr[,3][as.numeric(design.df$Trt)]
trt.contr4 = trt.contr[,4][as.numeric(design.df$Trt)]
trt.contr5 = trt.contr[,5][as.numeric(design.df$Trt)]


summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt",trt.contr = list(Tag = NA,
                                       Trt = list(  "1" = trt.contr1,
                                                    "2" = trt.contr2,
                                                    "3" = trt.contr3,
                                                    "4" = trt.contr4,
                                                    "5" = trt.contr5)))



trt.contr = test(X.trt = X.trt[bestInd,],  PP, Rep=trt.rep)$e.vec

trt.contr1 = trt.contr[,1][as.numeric(design.df$Trt)]
trt.contr2 = trt.contr[,2][as.numeric(design.df$Trt)]
trt.contr3 = trt.contr[,3][as.numeric(design.df$Trt)]
trt.contr4 = trt.contr[,4][as.numeric(design.df$Trt)]
trt.contr5 = trt.contr[,5][as.numeric(design.df$Trt)]
trt.contr6 = trt.contr[,6][as.numeric(design.df$Trt)]
trt.contr7 = trt.contr[,7][as.numeric(design.df$Trt)]


summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt",trt.contr = list(Tag = NA,
                                       Trt = list(  "1" = trt.contr1,
                                                    "2" = trt.contr2,
                                                    "3" = trt.contr3,
                                                    "4" = trt.contr4,
                                                    "5" = trt.contr5,
                                                    "6" = trt.contr6,
                                                    "7" = trt.contr7)))


info.mat =  matMulti( Pb - mK(n) , Z1.mat[bestInd,], C.cage)

if(any(abs(info.mat) > 1e-7)){
  PP = Pb - mK(n) - matMulti1(Pb - mK(n), ginv( matMulti(Pb - mK(n),
      Z1.mat[bestInd,], C.cage)), C.cage, Z1.mat[bestInd,])

  PP1 = matMulti1(PP, ginv( matMulti(PP , Z1.mat[bestInd,], C.ani)),
      C.ani, Z1.mat[bestInd,])

} else{
  PP1 = matMulti1(Pb - mK(n), ginv( matMulti((mI(n) - Pb),
      Z1.mat[bestInd,], C.ani)), C.ani, Z1.mat[bestInd,])
}

