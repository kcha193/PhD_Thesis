
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



test =
function(X.trt, blk.proj, C.trt.mat = diag(ncol(X.trt)), Rep){
  if( length(dim(X.trt))==3) X.trt = X.trt[,,1]

  info.mat = C.trt.mat %*% t(X.trt) %*% blk.proj %*% X.trt %*% C.trt.mat
  trace = tr(info.mat)
  e.va = Re(eigen(info.mat)$va)
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


newLETTERS = sort(sapply(strsplit(levels(interaction(LETTERS, LETTERS)), "\\."),
                function(x) paste(x, collapse = "", sep = "")))
#newLETTERS=LETTERS
sourceDir(path = "C:/Users/kcha193/My Dropbox/R functions for two-phase experiments")

#C++ matrix multiplication
code <-
  'arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat Z = Rcpp::as<arma::mat>(Z_);

  return wrap(Z*trans(Y) * X * Y * Z);'

matMulti <- cxxfunction( signature(X_="numeric",
Y_="numeric", Z_="numeric"), code, plugin="RcppArmadillo")

code <-
  'arma::mat W = Rcpp::as<arma::mat>(W_);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat Z = Rcpp::as<arma::mat>(Z_);

  return wrap(W*Z*Y * X * Y * trans(Z)*W);'

matMulti1 <- cxxfunction( signature(W_="numeric",X_="numeric",
Y_="numeric", Z_="numeric"), code, plugin="RcppArmadillo")

##############################################################################################
optRBD = function(nTrt, bRep, nCag, tRep, nPlot) {
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt

    nAni = nTrt * bRep

    if(!is.wholenumber( n/nPlot)){
     stop("The samples from Phase 1 experiment can not be fitted to the Phase 2 experiment!")
    }


    cat("Design parameters:\n")
    cat("Trt:", nTrt, "Ani:", nAni, "Tag:", nPlot, "Run:",
        nBlk, "\n")

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

print(summary.aov.onePhase(phase1DesignEX1, blk.str = "Cag/Ani", trt.str = "Trt"))

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

    fillIn = function(run, tag, phase1.mat, count, ani.limit, trt.limit) {

        for (i in 1:length(count)) {

            check1 = sum(phase1.mat[run, , 1] == names(count)[i]) < ani.limit[2]
            check2 = sum(phase1.mat[, tag, 1] == names(count)[i]) < ani.limit[1]
            trt = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == names(count)[i])])
            check3 = sum(phase1.mat[run, , 2] == trt) < trt.limit[2]
            check4 = sum(phase1.mat[, tag, 2] == trt) < trt.limit[1]

            if (length(count) == 1) {
                return(names(count)[i])
            }

            if (check1 && check2 && check3 && check4) {
                return(names(count)[i])
            }

        }
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
                phase1.mat[nBlk, j, 1] = fillIn(nBlk, j, phase1.mat, count, rep(Z1.rep,
                  2), c(nBlk/nTrt, trt.rep))

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

                  phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, rep(Z1.rep,
                    2), c(nBlk/nTrt, nTrt))

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
                  phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, ani.limit,
                    c(nBlk/nTrt, trt.rep/tRep))

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
    #print(phase1.mat)


    Z1.des = as.numeric(as.factor(t(phase1.mat[, , 1])))
    Z1.mat = matrix(0, ncol = nZ1, nrow = nZ1 * Z1.rep)
    Z1.mat[cbind(1:(nZ1 * Z1.rep), Z1.des)] = 1


    # Parameter's of treatment
    trt.rep = n/nTrt

    trt.des = as.numeric(phase1DesignEX1$Trt[match(newLETTERS[Z1.des], phase1DesignEX1$Ani)])
    X.trt = matrix(0, ncol = nTrt, nrow = nTrt * trt.rep)
    X.trt[cbind(1:(nTrt * trt.rep), trt.des)] = 1
    #############################################################################################

    getPair = function(x, n) {
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


  swap.new <- function(ind) {

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

    swap.stage1.new <- function(ind) {

        idx <- seq(2, length(ind) - 2)
        changepoints <- sample(idx, size = 2, replace = FALSE)
        # Two checks for omitting the swap of two indtical animals and animals in the
        # same block
        check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]],
            ])
        check2 = diff((changepoints - 1)%/%nPlot) != 0

        while (check1 || check2) {
            changepoints <- sample(idx, size = 2, replace = FALSE)
            check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]],
                ])
            check2 = diff((changepoints - 1)%/%nPlot) != 0

        }

        tmp <- ind[getPair(changepoints[1], length(ind) - 1)]
        ind[getPair(changepoints[1], length(ind) - 1)] <- ind[getPair(changepoints[2],
            length(ind) - 1)]
        ind[getPair(changepoints[2], length(ind) - 1)] <- tmp
        return(ind)
    }

    swap.stage2.new <- function(ind) {

        idx <- seq(2, length(ind) - 2)
        changepoints <- sample(idx, size = 2, replace = FALSE)
        # Two checks for omitting the swap of two indtical animals and animals in the
        # same block
        check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]],
            ])
        check2 = diff(changepoints%%nPlot) != 0

        while (check1 || check2) {
            changepoints <- sample(idx, size = 2, replace = FALSE)
            check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]],
                ])
            check2 = diff(changepoints%%nPlot) != 0
        }
        tmp <- ind[getPair(changepoints[1], length(ind) - 1)]
        ind[getPair(changepoints[1], length(ind) - 1)] <- ind[getPair(changepoints[2],
            length(ind) - 1)]
        ind[getPair(changepoints[2], length(ind) - 1)] <- tmp
        return(ind)
    }


################################################################################
#The objective function

resDF = tr(projMat(Z1.mat)) - (nCag - 1) - (nTrt -1) - 1
C.cage  = (mI(nCag) - mK(nCag)) %x%  mK(nAni/nCag) 
C.ani =   mI(nCag)%x%  (mI(nAni/nCag) - mK(nAni/nCag))
cage.Rep = n/nCag

C.ani =   mI(nCag)%x%  (mI(nAni/nCag) - mK(nAni/nCag))
ani.Rep = n/nAni
trt.Rep = trt.rep

try.obj.fun.new = function(ind) {

    #inforamtion in within runs and tags stratum while the cage information is eliminated
    info.mat = matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)
    PP = matMulti1(blk.proj, ginv(info.mat), C.cage,  Z1.mat[ind[-length(ind)],])
    
   # aveEff1 = tr(info.mat/cage.Rep)/tr(PP)
    
    #if(!isTRUE(all.equal(aveEff1, 1))) return(0)
     
    if(any(abs(info.mat) > 1e-7)){     #check for any information in the cage stratum
      PP = blk.proj - PP 
    } else{       #if there is no information in the cage stratum then no decomposition is required. If decomposition is performed here, 
                  # some very werid thing is going to happen. This is because taking the inverse of a very small number will become huge. 
                  # Hence, large information of the within block within runs can be extracted. 
      PP = blk.proj
    }
    
    #Maximise the amount of animal information in within runs and tags
    info.mat = matMulti(PP , Z1.mat[ind[-length(ind)],], C.ani)
    #e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    #can.eff = e.va[-which(e.va < 1e-07)]/ani.Rep 

    #aveEff1 = 100/mean(1/can.eff)
    
      
   #Maximise the amount of treatment information in within runs and tags

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
    
    
    aveEff =    ((tr(PP) - len)/resDF)  * ((resDF-1)/resDF)  +  (1/(mean(1/can.eff))) * (1/resDF)
    
    
    return(aveEff*100)
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

initialTemp = function(iter, newInit, swap){
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  U = numeric(iter)
  temp=matrix(0, nrow = iter, ncol = length(newInit))
  temp[1,] = newInit
  U[1] =   obj.fun.new(temp[1,])
  for(i in 2:iter){
    setTxtProgressBar(pb, i)
    temp[i,] = swap(temp[i-1,])
    U[i] = obj.fun.new(temp[i,])
  }
  close(pb)
  
  #if(any(U == 0)) U = U[-which(U ==0)]
  
  #y2 = max(U) - min(U) 
  
  #y2low = min(range((abs(diff(U))[which(abs(diff(U))>1e-5)])))
  if(all(U==0) || all(diff(U)<1e-14)){
    return(list(y2 = 100, y2low = 1e-5, newInit = newInit))
  
  }else{
    xxx = max(U, na.rm = TRUE) - range(U[which(U < (max(U, na.rm = TRUE)-1e-10))], na.rm = TRUE)
    
    if(isTRUE(all.equal(xxx[1], xxx[2]))) xxx[2] = 1e-5
     
    newInit = temp[which(U == max(U))[1],]
    return(list(y2 = max(xxx), y2low = min(xxx), newInit =newInit))
  }
}   

################################################################################
#Two-stages swapping method
twoStageRBD = function(init, iter = 50000) {
    print(c(nTrt, bRep, nCag, tRep, nPlot))
   
    newInit = c(init, init[1])

    old = obj.fun.new(newInit)

    cat("Level: 1, Finding the temperatures for both stages\n")

    temp = initialTemp(iter, newInit, swap.stage2.new)

    res.stage2 = temp$newInit
    y2.stage2 = temp$y2
    y2low.stage2 = temp$y2low

    cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")

    temp = initialTemp(iter, newInit, swap.stage1.new)

    res.stage1 = temp$newInit
    y2.stage1 = temp$y2
    y2low.stage1 = temp$y2low

    cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")

    temp = initialTemp(iter, newInit, swap.new)

    res = temp$newInit
    y2 = temp$y2
    y2low = temp$y2low
    
    cat("Stage 0 temperature range: ", y2, " to ", y2low, "\n")

    test.newInit = apply(rbind(res,res.stage1,res.stage2),1, obj.fun.new)
    newInit = rbind(res,res.stage1,res.stage2)[which(test.newInit == max(test.newInit))[1],]

    #newInit = if(obj.fun.new(res.stage1) > obj.fun.new(res.stage2)){ res.stage1} else{ res.stage2}

    cat("Level: 2, Current temp: ", y2.stage2, ",", y2.stage1, ",", y2, "\n")

    res <- optim(newInit, obj.fun.new, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
        tmax = 1000, trace = TRUE, REPORT = 100, fnscale = -1))

    res1 <- optim(res$par, obj.fun.new, swap.stage1.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
        tmax = 1000, trace = TRUE, REPORT = 100, fnscale = -1))

    res1 <- optim(res1$par, obj.fun.new, swap.new, method = "SANN", control = list(maxit = iter, temp = y2,
        tmax = 1000, trace = TRUE, REPORT = 100, fnscale = -1))

    cur = obj.fun.new(res1$par)

    print(test.obj.fun.new(res1$par))
    
    if(cur == 100) return(res1$par[-length(res1$par)])

    while ((cur - old) > 1e-06) {
        res1 <- optim(res1$par, obj.fun.new, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = 1000, trace = TRUE, REPORT = 100, fnscale = -1))

        res1 <- optim(res1$par, obj.fun.new, swap.stage1.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = 1000, trace = TRUE, REPORT = 100, fnscale = -1))

        res1 <- optim(res1$par, obj.fun.new, swap.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = 1000, trace = TRUE, REPORT = 100, fnscale = -1))

        old = cur
        cur = obj.fun.new(res1$par)

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
        
        old = obj.fun.new(res1$par)

        cat("Level: ", i + 2, ", Current temp: ", y2.stage2, ",", y2.stage1, ",", y2, "\n")
        res1 <- optim(res1$par, obj.fun.new, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = 1000, trace = TRUE, REPORT = 100, fnscale = -1))

        res1 <- optim(res1$par, obj.fun.new, swap.stage1.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = 1000, trace = TRUE, REPORT = 100, fnscale = -1))
        
        res1 <- optim(res1$par, obj.fun.new, swap.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
            tmax = 1000, trace = TRUE, REPORT = 100, fnscale = -1))

        cur = obj.fun.new(res1$par)
          print(test.obj.fun.new(res1$par))

        if(cur == 100) return(res1$par[-length(res1$par)])

        while ((cur - old) > 1e-06) {
            unstable = TRUE
            res1 <- optim(res1$par, obj.fun.new, swap.stage2.new, method = "SANN", control = list(maxit = iter,
                temp = y2.stage2, tmax = 1000, trace = TRUE, REPORT = 100, fnscale = -1))

            res1 <- optim(res1$par, obj.fun.new, swap.stage1.new, method = "SANN", control = list(maxit = iter,
                temp = y2.stage1, tmax = 1000, trace = TRUE, REPORT = 100, fnscale = -1))

            res1 <- optim(res1$par, obj.fun.new, swap.new, method = "SANN", control = list(maxit = iter, temp = y2.stage1,
              tmax = 1000, trace = TRUE, REPORT = 100, fnscale = -1))

            old = cur
            cur = obj.fun.new(res1$par)
               print(test.obj.fun.new(res1$par))

            if(cur == 100) return(res1$par[-length(res1$par)])

        }


        i = i + 1
    }

    print(c(nTrt, bRep, nCag, tRep, nPlot))

    return(res1$par[-length(res1$par)])
}

new.Z1.mat=  Z1.mat[twoStageRBD(1:n, iter = 1000),]

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

    cat("Phase 2 theoretical ANOVA:\n")
levels(design.df$Ani) = rep(1:(nlevels(design.df$Ani)/nCag), nCag)

print(summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Cag/Ani",
trt.str = "Tag + Trt") )

    return(design.df)
}


################################################################################
#design.summary

design.summary =
function(design.df, simple = TRUE) {


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

    C.cage = (mI(nCag) - mK(nCag)) %x% mK(nAni)
    C.ani = mI(nCag) %x% (mI(nAni) - mK(nAni))
    cage.Rep = n/nCag

    C.ani = mI(nCag) %x% (mI(nAni) - mK(nAni))

    Pb = projMat(Run.mat)
    Pb1 = projMat(Tag.mat)

    blk.proj = (mI(n) - Pb) %*% (mI(n) - Pb1)

    info.mat = matMulti(mI(n) - Pb, Ani.mat, C.cage)

    if (any(abs(info.mat) > 1e-07)) {
        PP = (mI(n) - Pb) - matMulti1((mI(n) - Pb), ginv(matMulti((mI(n) - Pb), Ani.mat,
            C.cage)), C.cage, Ani.mat)

        PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
    } else {
        PP1 = matMulti1((mI(n) - Pb), ginv(matMulti((mI(n) - Pb), Ani.mat, C.ani)),
            C.ani, Ani.mat)

    }


    cat("Design parameters:\n")
    cat("Trt:", ncol(Trt.mat), "Ani:", ncol(Ani.mat), "Cag:", ncol(Cag.mat), "Tag:", ncol(Tag.mat), "Run:",
        ncol(Run.mat), "\n")

    if(simple){
    cat("Animal design:\n")
    print(matrix(paste(design.df$Cag, design.df$Ani, sep = ""), nrow = nBlk, ncol = nPlot,
        byrow = TRUE))

    cat("Animal efficiency:\n")
    print(test(X.trt = mI(n), blk.proj = PP1, Rep = 1)[c("nCan", "can.eff")])

    info.mat = matMulti(blk.proj, Ani.mat, C.cage)

    if (any(abs(info.mat) > 1e-07)) {
        PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)),
            C.cage, Ani.mat)

        PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
    } else {
        PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani,
            Ani.mat)

    }


    cat("Treatment design:\n")
    print(matrix(design.df$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE))

    cat("Treatment incidence matrix:\n")
    print((N = with(design.df, table(Trt, Run))))

    cat("Treatment concurrence matrix:\n")
    print(N %*% t(N))

    cat("Treatment efficiency:\n")

    print(test(X.trt = Trt.mat, PP1, Rep = n/ncol(Trt.mat)))
    }

    cat("Phase 1 theoretical ANOVA:\n")
    print(summary.aov.onePhase(design.df, blk.str = "Cag/Ani", trt.str = "Trt"))

    cat("Phase 2 theoretical ANOVA:\n")
    print(summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Cag/Ani", trt.str = "Tag + Trt"))

}



design.df = optRBD(6,3,3,2,4)

design.df$Cag = as.factor(toupper(design.df$Cag))

design.summary.RBD(design.df, FALSE)


fileName = paste("RBDdesign_", 8, "trt_", 80, "ani_", 10,"blk_",  4, "tag_", 
40, "run.Rdata", sep = "")

fileName

save(design.df, file = fileName)
