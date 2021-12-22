#################################################################################
#Required functions
 
sourceDir(path = "C:/Users/kcha193/My Dropbox/R functions for two-phase experiments")

sourceDir <- function(path, trace = TRUE, ...) {
    library(MASS)
    library(lattice)
    library(inline)
    library(compiler)
    library(formatR)
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        if (trace)
            cat(nm, ":")
        source(file.path(path, nm), ...)
        if (trace)
            cat("\n")
    }
}



test.CRD = function(X.trt, blk.proj, C.trt.mat = diag(ncol(X.trt)), Rep) {
    if (length(dim(X.trt)) == 3)
        X.trt = X.trt[, , 1]

    info.mat = C.trt.mat %*% t(X.trt) %*% blk.proj %*% X.trt %*% C.trt.mat
    trace = tr(info.mat)
    e.va = eigen(info.mat)$va
    e.vec = eigen(info.mat)$vec
    can.eff = e.va[-which(e.va < 1e-07)]/Rep
    list(trace = trace, nCan = length(can.eff), can.eff = can.eff, ave.eff = 1/mean(1/can.eff),
        e.vec = e.vec)
}

multiLETTERS = function(x) {

    if (max(x) < 26) {
        LETTERS[x]
    } else {
        newLETTERS = sort(sapply(strsplit(levels(interaction(LETTERS, LETTERS)),
            "\\."), function(x) paste(x, collapse = "", sep = "")))
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

getPair = function(x, n, nPlot) {
    y = ceiling(x/nPlot)
    if (nPlot == 4) {
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
            } else if (y%%2 == 0) {
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

initialCRD = function(nTrt, bRep, tRep, nPlot) {
    newLETTERS = sort(sapply(strsplit(levels(interaction(LETTERS, LETTERS)), "\\."),
            function(x) paste(x, collapse = "", sep = "")))
            
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt
    
    nAni = nTrt * bRep
    
    if (!is.wholenumber(n/nPlot)) {
        stop("The samples from Phase 1 experiment can not be fitted to the Phase 2 experiment!")
    }
    
    
    # cat('Design parameters:\n') cat('Trt:', nTrt, 'Ani:', nAni, 'Tag:', nPlot, 'Run:', nBlk,
    # '\n')
    
    phase1DesignEX1 <- local({
        Ani = 1:nAni
        Trt = 1:nTrt
        
        data.frame(cbind(Ani, Trt))
    })
    
    phase1DesignEX1$Ani = as.factor(newLETTERS[phase1DesignEX1$Ani])
    phase1DesignEX1$Trt = as.factor(multiLetters(phase1DesignEX1$Trt))
    
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
    
    fillIn = function(run, tag, phase1.mat, count, ani.limit, trt.limit) {
        
        for (i in 1:length(count)) {
            
            check1 = sum(phase1.mat[run, , 1] == names(count)[i]) < ani.limit[2]
            check2 = sum(phase1.mat[, tag, 1] == names(count)[i]) < ani.limit[1]
            trt = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == names(count)[i])])
            check3 = sum(phase1.mat[run, , 2] == trt) < trt.limit[2]
            check4 = sum(phase1.mat[, tag, 2] == trt) < trt.limit[1]
            
            # print(c(check1, check2, check3, check4))
            
            if (length(count) == 1) {
                return(names(count)[i])
            }
            
            if (check1 && check2 && check3 && check4) {
                return(names(count)[i])
            } else if (run%%2 == 0 && tag%%2 == 1) {
                return(phase1.mat[run - 1, tag + 1, 1])
            } else if (run%%2 == 0 && tag%%2 == 0) {
                return(phase1.mat[run - 1, tag - 1, 1])
            } else if (run == nBlk) {
                return(names(count)[i])
            }
            
        }
        
        return(names(which(count == max(count)))[1])
    }
    
    ### fill-in then check ####
    ani.char = newLETTERS[1:nZ1]
    count = rep(Z1.rep, nZ1)
    names(count) = ani.char
    
    phase1.mat = array("-1", c(nBlk, nPlot, 2))
    
    len = (nBlk%/%tRep) * tRep
    
    ani.limit = c(len, nPlot)/nZ1
    trt.limit = c(len, nPlot)/nTrt
    
    for (i in 1:len) {
        for (j in 1:nPlot) {
            
            phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, ani.limit, trt.limit)
            
            phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == 
                phase1.mat[i, j, 1])])
            
            count[phase1.mat[i, j, 1]] = count[phase1.mat[i, j, 1]] - 1
            
            # count = sort(count, d=TRUE)
            if (any(count == 0)) 
                count = count[-which(count == 0)]
        }
    }
    
    if (nBlk != len) {
        if (nBlk%%tRep == 1) {
            for (j in 1:nPlot) {
                phase1.mat[nBlk, j, 1] = fillIn(nBlk, j, phase1.mat, count, rep(Z1.rep, 2), c(nBlk/nTrt, 
                  trt.rep))
                
                phase1.mat[nBlk, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == 
                  phase1.mat[nBlk, j, 1])])
                count[phase1.mat[nBlk, j, 1]] = count[phase1.mat[nBlk, j, 1]] - 1
                
                # count = sort(count, d=TRUE)
                if (any(count == 0)) 
                  count = count[-which(count == 0)]
            }
            
        } else if (nBlk%%tRep == 2) {
            for (i in (len + 1):nBlk) {
                for (j in 1:nPlot) {
                  
                  phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, rep(Z1.rep, 2), c(nBlk/nTrt, 
                    nTrt))
                  
                  phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == 
                    phase1.mat[i, j, 1])])
                  
                  count[phase1.mat[i, j, 1]] = count[phase1.mat[i, j, 1]] - 1
                  
                  # count = sort(count, d=TRUE)
                  if (any(count == 0)) 
                    count = count[-which(count == 0)]
                }
            }
        } else {
            for (i in (len + 1):nBlk) {
                for (j in 1:nPlot) {
                  
                  ani.limit = c(length((len + 1):nBlk), nPlot)/length(count)
                  phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, ani.limit, c(nBlk/nTrt, 
                    trt.rep/tRep))
                  
                  phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == 
                    phase1.mat[i, j, 1])])
                  
                  count[phase1.mat[i, j, 1]] = count[phase1.mat[i, j, 1]] - 1
                  
                  # count = sort(count, d=TRUE)
                  if (any(count == 0)) 
                    count = count[-which(count == 0)]
                }
            }
            
        }
    }
    # print(phase1.mat)
    
    as.numeric(as.factor(t(phase1.mat[, , 1])))
}

###################################################################################

   swap.stage1.new <- function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt) {
        
        idx <- seq(2, length(ind) - 2)
        changepoints <- sample(idx, size = 2, replace = FALSE)
        # Two checks for omitting the swap of two indtical animals and animals in the same block
        check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]], ])
        check2 = diff((changepoints - 1)%/%nPlot) != 0
        
        while (check1 || check2) {
            changepoints <- sample(idx, size = 2, replace = FALSE)
            check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]], ])
            check2 = diff((changepoints - 1)%/%nPlot) != 0
            
        }
        
        tmp <- ind[getPair(changepoints[1], length(ind) - 1, nPlot)]
        ind[getPair(changepoints[1], length(ind) - 1, nPlot)] <- ind[getPair(changepoints[2], length(ind) - 
            1, nPlot)]
        ind[getPair(changepoints[2], length(ind) - 1, nPlot)] <- tmp
        return(ind)
    }
    
    swap.stage2.new <- function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt) {
        
        idx <- seq(2, length(ind) - 2)
        changepoints <- sample(idx, size = 2, replace = FALSE)
        # Two checks for omitting the swap of two indtical animals and animals in the same block
        check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]], ])
        check2 = diff(changepoints%%nPlot) != 0
        
        while (check1 || check2) {
            changepoints <- sample(idx, size = 2, replace = FALSE)
            check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]], ])
            check2 = diff(changepoints%%nPlot) != 0
        }
        tmp <- ind[getPair(changepoints[1], length(ind) - 1, nPlot)]
        ind[getPair(changepoints[1], length(ind) - 1, nPlot)] <- ind[getPair(changepoints[2], length(ind) - 
            1, nPlot)]
        ind[getPair(changepoints[2], length(ind) - 1, nPlot)] <- tmp
        return(ind)
    }
    
    
    swap.stage3.new <- function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt) {
        
        idx <- seq(2, length(ind) - 2)
        changepoints <- sample(idx, size = 2, replace = FALSE)
        # Two checks for omitting the swap of two indtical animals and animals in the same block
        check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]], ])
        check2 = diff((changepoints - 1)%/%nPlot) == 0
        check3 = diff(changepoints%%nPlot) == 0
        
        while (check1 || check2 || check3) {
            changepoints <- sample(idx, size = 2, replace = FALSE)
            check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]], ])
            check2 = diff((changepoints - 1)%/%nPlot) == 0
            check3 = diff(changepoints%%nPlot) == 0
        }
        tmp <- ind[getPair(changepoints[1], length(ind) - 1, nPlot)]
        ind[getPair(changepoints[1], length(ind) - 1, nPlot)] <- ind[getPair(changepoints[2], length(ind) - 
            1, nPlot)]
        ind[getPair(changepoints[2], length(ind) - 1, nPlot)] <- tmp
        return(ind)
    }
    
    
    
    
    
    
    ######################################################################################
    try.obj.fun.CRD = function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt) {
        
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
        
        aveEff = (2/3) * aveEff + 1/3 * ((len + 1/mean(1/can.eff))/nTrt)
        
        return(aveEff)
    }
    
    obj.fun.CRD = function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt) {
        
        ans = try(try.obj.fun.CRD(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt), silent = TRUE)
        
        ifelse(class(ans) == "try-error", 0, ans)
    }
    
    
    ################################################################################ Initialise
    ################################################################################ the
    ################################################################################ temperature
    
    test.obj.fun.CRD = function(ind, Z1.mat, X.trt, Z1.rep, trt.rep, blk.proj) {
        
        
        # inforamtion in within runs and tags stratum while the cage information is eliminated
        
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
######################################################################################## 

initialTemp = function(iter, newInit, swap, obj.fun, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt){
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  U = numeric(iter)
  temp = matrix(0, nrow = iter, ncol = length(newInit))
  temp[1,] = newInit
  U[1] =   obj.fun(temp[1,], Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)
  tol = 10
  check.U = numeric(tol) 
  
  for(i in 2:iter){
    setTxtProgressBar(pb, i)
    
    temp[i,] = swap(temp[i-1,], Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)
    U[i] = obj.fun(temp[i,], Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)
    
    #check.U = c(check.U, U[i])
    #check.U = check.U[-1]
   
    if(all(outer(check.U,check.U,"=="))) temp[i,] = newInit
  }
  close(pb)
  
  if(all(U==0) || all(diff(U)<1e-8)){
    return(list(y2 = 100, y2low = 1e-5, newInit = newInit))
  
  }else{
    xxx = max(U, na.rm = TRUE) - range(U[which(U < (max(U, na.rm = TRUE)-1e-10))], na.rm = TRUE)
    
    if(isTRUE(all.equal(xxx[1], xxx[2]))) xxx[2] = 1e-5
     
    newInit = temp[which(U == max(U))[1],]
    return(list(y2 = max(xxx), y2low = min(xxx), newInit =newInit))
  }
}   

################################################################################


optThreeStage = function(init, iter, obj.fun, test.obj.fun, swap.stage1.new, swap.stage2.new, 
    swap.stage3.new, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt) {
    
   newInit = c(init, init[1])
    #browser()
    old = obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

    cat("Level: 1, Finding the temperatures for both stages\n")

    temp = initialTemp(iter, newInit, swap.stage2.new, obj.fun, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

    res.stage2 = temp$newInit
    y2.stage2 = temp$y2
    y2low.stage2 = temp$y2low

    temp = initialTemp(iter, newInit, swap.stage1.new, obj.fun, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

    res.stage1 = temp$newInit
    y2.stage1 = temp$y2
    y2low.stage1 = temp$y2low

    temp = initialTemp(iter, newInit, swap.stage3.new, obj.fun, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

    res = temp$newInit
    y2 = temp$y2
    y2low = temp$y2low
    
    cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
    cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
    cat("Stage 0 temperature range: ", y2, " to ", y2low, "\n")
 
    test.newInit = apply(rbind(res,res.stage1,res.stage2),1, obj.fun, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)
    newInit = rbind(res,res.stage1,res.stage2)[which(test.newInit == max(test.newInit))[1],]
    
    cur = obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)
    print(test.obj.fun(newInit, Z1.mat, X.trt, Z1.rep, trt.rep, blk.proj))
    
    while ((cur - old) > 1e-06) {
        
        temp = initialTemp(iter, newInit, swap.stage2.new, obj.fun, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)
    
        res.stage2 = temp$newInit
        y2.stage2 = temp$y2
        y2low.stage2 = temp$y2low
    
    
        temp = initialTemp(iter, newInit, swap.stage1.new, obj.fun, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)
    
        res.stage1 = temp$newInit
        y2.stage1 = temp$y2
        y2low.stage1 = temp$y2low
    
    
        temp = initialTemp(iter, newInit, swap.stage3.new, obj.fun, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)
    
        res = temp$newInit
        y2 = temp$y2
        y2low = temp$y2low
           
        test.newInit = apply(rbind(res,res.stage1,res.stage2),1, obj.fun, Z1.mat, X.trt, Z1.rep, trt.rep, blk.proj, nTrt)
        newInit = rbind(res,res.stage1,res.stage2)[which(test.newInit == max(test.newInit))[1],]

        old = cur
        
        cur = obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)
        print(test.obj.fun.new(newInit, Z1.mat, X.trt, Z1.rep, trt.rep, blk.proj))     
   
        cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
        cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
        cat("Stage 0 temperature range: ", y2, " to ", y2low, "\n")
    
     }
      

    old = obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)
  
    cat("Level: 2, Current temp: ", y2.stage2, ",", y2.stage1, ",", y2, "\n")

    res <- optim(newInit, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
        tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

    res1 <- optim(res$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
        tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

    res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, temp = y2,
        tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

    cur = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

    print(test.obj.fun(res1$par, Z1.mat, X.trt, Z1.rep, trt.rep, blk.proj))
    
    if(cur == 1) return(res1$par[-length(res1$par)])

    while ((cur - old) > 1e-06) {
        res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

        res1 <- optim(res1$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

        res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, temp = y2,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

        old = cur
        cur = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

        print(test.obj.fun(res1$par, Z1.mat, X.trt, Z1.rep, trt.rep, blk.proj))


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
        
        old = obj.fun(res1$par, Z1.mat, X.trt, Z1.rep, trt.rep, blk.proj)

        cat("Level: ", i + 2, ", Current temp: ", y2.stage2, ",", y2.stage1, ",", y2, "\n")
        res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

        res1 <- optim(res1$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)
        
        res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, temp = y2,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

        cur = obj.fun(res1$par, Z1.mat, X.trt, Z1.rep, trt.rep, blk.proj, nTrt)
          print(test.obj.fun(res1$par, Z1.mat, X.trt, Z1.rep, trt.rep, blk.proj))

        if(cur == 1) return(res1$par[-length(res1$par)])

        while ((cur - old) > 1e-06) {
            unstable = TRUE
            res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter,
                temp = y2.stage2, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
                Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

            res1 <- optim(res1$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter,
                temp = y2.stage1, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
                Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

            res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, temp = y2,
              tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
              Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)

            old = cur
            cur = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt)
            print(test.obj.fun(res1$par, Z1.mat, X.trt, Z1.rep, trt.rep, blk.proj))

            if(cur == 1) return(res1$par[-length(res1$par)])

        }


        i = i + 1
    }

 
    return(res1$par[-length(res1$par)])
}

################################################################################

optCRD = function(nTrt, bRep, tRep, nPlot, iter = 10000) {
    newLETTERS = sort(sapply(strsplit(levels(interaction(LETTERS, LETTERS)), "\\."),
      function(x) paste(x, collapse = "", sep = "")))
    
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
        Ani = 1:nAni
        Trt = 1:nTrt
        
        data.frame(cbind(Ani, Trt))
    })
    phase1DesignEX1$Ani = as.factor(newLETTERS[phase1DesignEX1$Ani])
    phase1DesignEX1$Trt = as.factor(multiLetters(phase1DesignEX1$Trt))
    
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
    
    Z1.des = initialCRD(nTrt = nTrt, bRep = bRep, tRep = tRep, nPlot = nPlot)
    
    Z1.mat = matrix(0, ncol = nZ1, nrow = nZ1 * Z1.rep)
    Z1.mat[cbind(1:(nZ1 * Z1.rep), Z1.des)] = 1
    
    
    # Parameter's of treatment
    trt.rep = n/nTrt
    
    trt.des = as.numeric(phase1DesignEX1$Trt[match(newLETTERS[Z1.des], phase1DesignEX1$Ani)])
    X.trt = matrix(0, ncol = nTrt, nrow = nTrt * trt.rep)
    X.trt[cbind(1:(nTrt * trt.rep), trt.des)] = 1
   
     
    ################################################################################ Initialise
    ################################################################################ the
    ################################################################################ temperature
    
    new.Z1.mat = Z1.mat[optThreeStage(init = 1:n, iter = iter, obj.fun = obj.fun.CRD, 
        test.obj.fun = test.obj.fun.CRD, swap.stage1.new = swap.stage1.new, 
        swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new, 
        Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, trt.rep = trt.rep,
        blk.proj = blk.proj, nTrt = nTrt), 
        ]
    
    ################################################################################ Set the
    ################################################################################ design from
    ################################################################################ the search
    
    colnames(new.Z1.mat) = sort(levels(interaction(newLETTERS[1:nZ1])))
    new.Z1.des = apply(new.Z1.mat, 1, function(x) colnames(new.Z1.mat)[which(as.logical(x))])
    new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))
    
    design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))
    
    design.df = cbind(design.df, Ani = t(new.Z1.des), Trt = phase1DesignEX1[match(as.character(t(new.Z1.des)), 
        as.character(phase1DesignEX1$Ani)), ]$Trt)
     
    return(design.df)
} 

################################################################################

design.df = optCRD(6,3,2,4, 1000)


################################################################################
#design summary for CRD

design.summary.CRD = function(design.df, simple = TRUE) {

    n = nrow(design.df)
    nBlk = nlevels(design.df$Run)
    nPlot = nlevels(design.df$Tag)
    Run.mat = with(design.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
    Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
    Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))

    Pb = projMat(Run.mat)
    Pb1 = projMat(Tag.mat)

    blk.proj = (mI(n) - Pb)

    cat("Design parameters:\n")
    cat("Trt:", ncol(Trt.mat), "Ani:", ncol(Ani.mat), "Tag:", ncol(Tag.mat), "Run:",
        ncol(Run.mat), "\n")

    if(simple){
    cat("Animal design:\n")
    print(matrix(design.df$Ani, nrow = nBlk, ncol = nPlot, byrow = TRUE))

    cat("Animal incidence matrix:\n")
    print((N = with(design.df, table(Ani, Run))))

    cat("Animal concurrence matrix:\n")
    print(N %*% t(N))

    cat("Animal efficiency:\n")
    print(test.CRD(X.trt = Ani.mat, blk.proj, Rep = n/ncol(Ani.mat)))

    cat("Treatment design:\n")
    print(matrix(design.df$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE))

    cat("Treatment incidence matrix:\n")
    print((N = with(design.df, table(Trt, Run))))

    cat("Treatment concurrence matrix:\n")
    print(N %*% t(N))

    cat("Treatment efficiency:\n")
    print(test.CRD(X.trt = Trt.mat, (mI(n) - Pb) %*% (mI(n) - Pb1), Rep = n/ncol(Trt.mat)))
    }

    cat("Phase 1 theoretical ANOVA:\n")
    print(summary.aov.onePhase(design.df, blk.str = "Ani", trt.str = "Trt"))

    cat("Phase 2 theoretical ANOVA:\n")
    print(summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt"))


}





design.summary.CRD(design.df, FALSE)

design.summary.CRD(design.df)
