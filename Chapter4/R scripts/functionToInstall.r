
matMulti <- function(P, X, C) C %*% t(X) %*% P %*% X %*% C

####################################################################################

matMulti1 <- function(P, G, C, X)  P %*% X %*% C %*% G %*% C %*% t(X) %*% P 

###################################################################################

initialRBD.nolimit = function(nTrt, bRep, nCag, tRep, nPlot) {
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt

    nAni = nTrt * bRep

    if (!is.wholenumber(n/nPlot)) {
        stop("The samples from Phase 1 experiment can not be fitted to the Phase 2 experiment!")
    }


    phase1DesignEX1 <- local({

        Cag = rep(1:nCag, each = (nTrt * bRep)/nCag)
        Ani = 1:nAni
        Trt = 1:nTrt

        data.frame(cbind(Cag, Ani, Trt))
    })

   phase1DesignEX1$Ani = as.factor(multiLetters(phase1DesignEX1$Ani))
    phase1DesignEX1$Trt = as.factor(multiLetters(phase1DesignEX1$Trt, FALSE))
    phase1DesignEX1$Cag = as.factor(phase1DesignEX1$Cag)



    nBlk = n/nPlot

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
    ani.char = multiLetters(1:nZ1)
    count = rep(Z1.rep, nZ1)
    names(count) = ani.char

    phase1.mat = array("-1", c(nBlk, nPlot, 2))

    len = (nBlk%/%tRep) * tRep

    ani.limit = c(len, nPlot)/nZ1
    trt.limit =c(Inf, Inf) # c(len, nPlot)/nTrt   

    for (i in 1:len) {
        for (j in 1:nPlot) {

            phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, ani.limit, trt.limit)

            phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[i,
                j, 1])])

            count[phase1.mat[i, j, 1]] = count[phase1.mat[i, j, 1]] - 1

            # count = sort(count, d=TRUE)
            if (any(count == 0))
                count = count[-which(count == 0)]
        }
    }

    if (nBlk != len) {
        if (nBlk%%tRep == 1) {
            for (j in 1:nPlot) {
                phase1.mat[nBlk, j, 1] = fillIn(nBlk, j, phase1.mat, count, rep(Z1.rep, 2), c(nBlk/nTrt, trt.rep))

                phase1.mat[nBlk, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[nBlk,
                  j, 1])])
                count[phase1.mat[nBlk, j, 1]] = count[phase1.mat[nBlk, j, 1]] - 1

                # count = sort(count, d=TRUE)
                if (any(count == 0))
                  count = count[-which(count == 0)]
            }

        } else if (nBlk%%tRep == 2) {
            for (i in (len + 1):nBlk) {
                for (j in 1:nPlot) {

                  phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, rep(Z1.rep, 2), c(nBlk/nTrt, nTrt))

                  phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[i,
                    j, 1])])

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
                  phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, ani.limit, c(nBlk/nTrt, trt.rep/tRep))

                  phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[i,
                    j, 1])])

                  count[phase1.mat[i, j, 1]] = count[phase1.mat[i, j, 1]] - 1

                  # count = sort(count, d=TRUE)
                  if (any(count == 0))
                    count = count[-which(count == 0)]
                }
            }

        }
    }

    return(as.numeric(as.factor(t(phase1.mat[,,1]))))
}

####################################################################################

orderCage = function(nTrt, bRep, nCag, tRep, nPlot, nAni) {
    test.Z1.des = initialRBD.nolimit(nTrt = nTrt, bRep = bRep, nCag = nCag, tRep = tRep, nPlot = nPlot)

    test.nBlk = length(test.Z1.des)/nPlot

    nAni = max(test.Z1.des)
    reSort.Z1.des = factor(test.Z1.des)
    temp = character()

    inc = nPlot/nCag
    start1 = seq(1, nAni/nCag, by = inc)


    for (i in 1:length(start1)){
        start2 = seq(start1[i], nAni, by = nAni/nCag)

        if ((test.nBlk%%2 == 1) && (i == length(start1))) {

            temp = c(temp, as.character((sapply(start2, function(x) seq(x, x + (inc/2) - 1, by = 1)))))
        } else {
            temp = c(temp, as.character(sapply(start2, function(x) seq(x, x + inc - 1, by = 1))))
        }
    }

    levels(reSort.Z1.des) = temp

    return(reSort.Z1.des)
}

# Allocation the cages Tag 4 Check for clostest tag 2, 4, and 8 Final function
sortCage =  function(nTrt, bRep, nCag, tRep, nPlot, nAni) {
if (nPlot == 8) {

    check1 = nCag - c(2, 4, 8)

    test.nCag = c(2, 4, 8)[which((nCag - c(2, 4, 8)) == min(check1[check1 > -0.01]))]

    dif = nCag - test.nCag
    test.bRep = bRep - (bRep/nCag * dif)
    test.nAni = nTrt * test.bRep


    reSort.Z1.des = orderCage(nTrt = nTrt, bRep = test.bRep, nCag = test.nCag, tRep = tRep, nPlot = nPlot, nAni = test.nAni)


    if (dif == 1) {
        Z1.des = c(as.numeric(as.character(reSort.Z1.des)), Z1.des[-(1:length(reSort.Z1.des))])

    } else if (dif == 2) {
        test2.bRep = (nAni - test.nAni)/nTrt

        reSort1.Z1.des = orderCage(nTrt = nTrt, bRep = test2.bRep, nCag = test2.bRep, tRep = tRep, nPlot = nPlot, nAni = nAni -
            test.nAni)

        Z1.des = c(as.numeric(as.character(reSort.Z1.des)), max(as.numeric(as.character(reSort.Z1.des))) + as.numeric(as.character(reSort1.Z1.des)))

    } else if (dif == 3) {
        test2.bRep = (nAni - test.nAni)/nTrt - 1

        reSort1.Z1.des = orderCage(nTrt = nTrt, bRep = test2.bRep, nCag = test2.bRep, tRep = tRep, nPlot = nPlot, nAni = nAni -
            test.nAni)


        Z1.des = c(as.numeric(as.character(reSort.Z1.des)), max(as.numeric(as.character(reSort.Z1.des))) + as.numeric(as.character(reSort1.Z1.des)),
            Z1.des[-(1:length(c(as.numeric(as.character(reSort.Z1.des)), as.numeric(as.character(reSort1.Z1.des)))))])

    } else {
        Z1.des = as.numeric(as.character(reSort.Z1.des))
    }



} else if (nPlot == 4) {
    #bRep = 8
    
    check1 = nCag - c(2, 4)

    test.nCag = c(2, 4)[which((nCag - c(2, 4)) == min(check1[check1 > -0.01]))]

    dif = nCag - test.nCag
    test.bRep = bRep - (bRep/nCag * dif)
    test.nAni = nTrt * test.bRep


    reSort.Z1.des = orderCage(nTrt = nTrt, bRep = test.bRep, nCag = test.nCag, tRep = tRep, nPlot = nPlot, nAni = test.nAni)


    if (dif == 1) {
        Z1.des = c(as.numeric(as.character(reSort.Z1.des)), Z1.des[-(1:length(reSort.Z1.des))])

    } else if (dif == 2 || dif == 4) {
        test2.bRep = (nAni - test.nAni)/nTrt

        reSort1.Z1.des = orderCage(nTrt = nTrt, bRep = test2.bRep, nCag = test2.bRep, tRep = tRep, nPlot = nPlot, nAni = nAni -
            test.nAni)

        Z1.des = c(as.numeric(as.character(reSort.Z1.des)), max(as.numeric(as.character(reSort.Z1.des))) + as.numeric(as.character(reSort1.Z1.des)))

    } else if (dif == 3 || dif == 5) {
        test2.bRep = (nAni - test.nAni)/nTrt - 1

        reSort1.Z1.des = orderCage(nTrt = nTrt, bRep = test2.bRep, nCag = test2.bRep, tRep = tRep, nPlot = nPlot, nAni = nAni -
            test.nAni)


        Z1.des = c(as.numeric(as.character(reSort.Z1.des)), max(as.numeric(as.character(reSort.Z1.des))) + as.numeric(as.character(reSort1.Z1.des)),
            Z1.des[-(1:length(c(as.numeric(as.character(reSort.Z1.des)), as.numeric(as.character(reSort1.Z1.des)))))])

    } else if (dif == 6) {
        test2.bRep = (nAni - test.nAni)/nTrt - 2

        reSort1.Z1.des = orderCage(nTrt = nTrt, bRep = test2.bRep, nCag = test2.bRep, tRep = tRep, nPlot = nPlot, nAni = nAni -
            test.nAni)

        test3.bRep = 2

        reSort2.Z1.des = orderCage(nTrt = nTrt, bRep = test3.bRep, nCag = test3.bRep, tRep = tRep, nPlot = nPlot, nAni = nAni -
            test.nAni)

        Z1.des = c(as.numeric(as.character(reSort.Z1.des)), max(as.numeric(as.character(reSort.Z1.des))) + as.numeric(as.character(reSort1.Z1.des)),
            max(as.numeric(as.character(reSort.Z1.des))) + max(as.numeric(as.character(reSort1.Z1.des))) + as.numeric(as.character(reSort2.Z1.des)))

    } else {
        Z1.des = as.numeric(as.character(reSort.Z1.des))
    }

}
 return(Z1.des)
}


###################################################################################

initialRBD = function(nTrt, bRep, nCag, tRep, nPlot) {
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt

    nAni = nTrt * bRep

    if (!is.wholenumber(n/nPlot)) {
        stop("The samples from Phase 1 experiment can not be fitted to the Phase 2 experiment!")
    }


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

    nBlk = n/nPlot


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

            phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[i,
                j, 1])])

            count[phase1.mat[i, j, 1]] = count[phase1.mat[i, j, 1]] - 1

            # count = sort(count, d=TRUE)
            if (any(count == 0))
                count = count[-which(count == 0)]
        }
    }

    if (nBlk != len) {
        if (nBlk%%tRep == 1) {
            for (j in 1:nPlot) {
                phase1.mat[nBlk, j, 1] = fillIn(nBlk, j, phase1.mat, count, rep(Z1.rep, 2), c(nBlk/nTrt, trt.rep))

                phase1.mat[nBlk, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[nBlk,
                  j, 1])])
                count[phase1.mat[nBlk, j, 1]] = count[phase1.mat[nBlk, j, 1]] - 1

                # count = sort(count, d=TRUE)
                if (any(count == 0))
                  count = count[-which(count == 0)]
            }

        } else if (nBlk%%tRep == 2) {
            for (i in (len + 1):nBlk) {
                for (j in 1:nPlot) {

                  phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, rep(Z1.rep, 2), c(nBlk/nTrt, nTrt))

                  phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[i,
                    j, 1])])

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
                  phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, ani.limit, c(nBlk/nTrt, trt.rep/tRep))

                  phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[i,
                    j, 1])])

                  count[phase1.mat[i, j, 1]] = count[phase1.mat[i, j, 1]] - 1

                  # count = sort(count, d=TRUE)
                  if (any(count == 0))
                    count = count[-which(count == 0)]
                }
            }

        }
    }

    return(as.numeric(as.factor(t(phase1.mat[,,1]))))
}


#################################################################################################


optRBD = function(nTrt, bRep, nCag, tRep, nPlot, iter = 10000) {
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt
    
    nAni = nTrt * bRep
    
    if (!is.wholenumber(n/nPlot)) {
        stop("The samples from Phase 1 experiment can not be fitted to the Phase 2 experiment!")
    }
     
    cat("Design parameters:\n")
    cat("Trt:", nTrt, "Ani:", nAni, "Cag:", nCag, "Tag:", nPlot, "Run:", nBlk, "\n")
    
    phase1DesignEX1 <- local({
        
        Cag = rep(1:nCag, each = (nTrt * bRep)/nCag)
        Ani = 1:nAni
        Trt = 1:nTrt
        
        data.frame(cbind(Cag, Ani, Trt))
    })
    
    phase1DesignEX1$Ani = as.factor(multiLetters(phase1DesignEX1$Ani))
    phase1DesignEX1$Trt = as.factor(multiLetters(phase1DesignEX1$Trt, FALSE))
    phase1DesignEX1$Cag = as.factor(phase1DesignEX1$Cag)
    
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
    
    if((nPlot %% nCag == 0) && (nPlot != nCag )){
      Z1.des = as.numeric(sortCage(nTrt = nTrt, bRep = bRep, nCag = nCag, tRep = tRep, 
                        nPlot = nPlot, nAni = nAni))
    }else {
      Z1.des = as.numeric(initialRBD(nTrt = nTrt, bRep = bRep, nCag = nCag, tRep = tRep, 
                        nPlot = nPlot))
    }
    
    
    Z1.mat = matrix(0, ncol = nZ1, nrow = nZ1 * Z1.rep)
    Z1.mat[cbind(1:(nZ1 * Z1.rep), Z1.des)] = 1
    
    
    # Parameter's of treatment
    trt.rep = n/nTrt
    
    trt.des = as.numeric(phase1DesignEX1$Trt[match(multiLetters(Z1.des), phase1DesignEX1$Ani)])
    X.trt = matrix(0, ncol = nTrt, nrow = nTrt * trt.rep)
    X.trt[cbind(1:(nTrt * trt.rep), trt.des)] = 1
    #############################################################################################
    #Contrast matrix for cage and animals 
    C.cage  = (mI(nCag) - mK(nCag)) %x%  mK(nAni/nCag) 
    cage.Rep = n/nCag

    C.ani =   mI(nCag)%x%  (mI(nAni/nCag) - mK(nAni/nCag))
    ani.Rep = n/nAni

    print("Step 1 optimisation:")
    firstStep = optThreeStage(init = 1:n, iter = iter/10, obj.fun = obj.fun.RBD, 
        test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage1.new, 
        swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new, 
        Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
        trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
        cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = NA)
       
    newResDF = test.obj.fun.RBD(c(firstStep,1), Z1.mat, X.trt, nPlot, Z1.rep, 
                              trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, 
                              ani.Rep, newResDF) $res
    
    print("Step 2 optimisation:")
    secondStep = optThreeStage(init = firstStep, iter = iter, obj.fun = obj.fun.RBD1, 
        test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage1.new, 
        swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new, 
        Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
        trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
        cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = newResDF)
        
    ans <- readline("Further improve the design (Y/N)?:")    
    
    while(ans == "Y"){
        newResDF = newResDF -1
        secondStep.new = optThreeStage(init = secondStep, iter = iter, obj.fun = obj.fun.RBD1, 
              test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage1.new, 
              swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new, 
              Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
              trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
              cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = newResDF)
        
        ans <- readline("Further improve the design (Y/N)?:")
        
        
        if(ans == "N") { 
                                                   
         if(as.integer(test.obj.fun.RBD(c(secondStep.new,1), Z1.mat, X.trt, nPlot, Z1.rep, 
                              trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, 
                              ani.Rep, newResDF)$aveEff) > 
            as.integer(test.obj.fun.RBD(c(secondStep,1), Z1.mat, X.trt, nPlot, Z1.rep, 
                              trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, 
                              ani.Rep, newResDF)$aveEff)){
                              
            secondStep  = secondStep.new
         } 
          
          break
        }
            
        secondStep = secondStep.new
    }
    
    new.Z1.mat = Z1.mat[secondStep,]
      
    
    ###############################################################################
    # Set the design from the search
    
    colnames(new.Z1.mat) = sort(levels(interaction(multiLetters(1:nZ1))))
    new.Z1.des = apply(new.Z1.mat, 1, function(x) colnames(new.Z1.mat)[which(as.logical(x))])
    new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))
    
    design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))
    
    design.df = cbind(design.df, 
              Cag = phase1DesignEX1[match(as.character(t(new.Z1.des)), 
                                    as.character(phase1DesignEX1$Ani)), ]$Cag, 
              Ani = t(new.Z1.des), 
              Trt = phase1DesignEX1[match(as.character(t(new.Z1.des)), 
                                    as.character(phase1DesignEX1$Ani)), ]$Trt)
    
    levels(design.df$Ani) = multiLetters(rep(1:(nlevels(design.df$Ani)/nCag), nCag) )

    return(design.df)
} 

#################################################################################
#First objective function

try.obj.fun.RBD = function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {

    #inforamtion in within runs and tags stratum while the cage information is eliminated
    info.mat = matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)
    PP = matMulti1(blk.proj, ginv(info.mat), C.cage,  Z1.mat[ind[-length(ind)],])
     
    if(any(abs(info.mat) > 1e-7)){     #check for any information in the cage stratum
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
    info.mat = matMulti(PP , Z1.mat[ind[-length(ind)],], C.ani)
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    can.eff = e.va[-which(e.va < 1e-07)]/ani.Rep 

    aveEff1 = 1/mean(1/can.eff)
     
    if(is.na(aveEff1)) aveEff1 = 0
   
    if(!isTRUE(all.equal(aveEff1, 1))) return(0)
     

    PP =  matMulti1(PP, ginv(info.mat), C.ani, Z1.mat[ind[-length(ind)],]) 

    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*%  PP %*% newX
    #len = length(can.eff)

    #len = qr(info.mat)$rank 
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
 
    can.eff = e.va[-which(e.va < 1e-07)]/trt.rep
    len = length(can.eff)
  
    if(!isTRUE(all.equal(len/(nTrt-1), 1))) return(0)
 
    aveEff =  tr(PP) - len 
    
    return(aveEff)
}


obj.fun.RBD = function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {

  ans = try(try.obj.fun.RBD(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
              nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF), silent = TRUE)
  
  ifelse(class(ans) == "try-error", 0, ans)
}

#####################################################################################
#Second objective function

try.obj.fun.RBD1 = function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {
       
    #inforamtion in within runs and tags stratum while the cage information is eliminated
    info.mat = matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)
    PP = matMulti1(blk.proj, ginv(info.mat), C.cage,  Z1.mat[ind[-length(ind)],])
     
    if(any(abs(info.mat) > 1e-7)){     #check for any information in the cage stratum
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
    info.mat = matMulti(PP , Z1.mat[ind[-length(ind)],], C.ani)
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    can.eff = e.va[-which(e.va < 1e-07)]/ani.Rep 

    aveEff1 = 1/mean(1/can.eff)
     
    if(is.na(aveEff1)) aveEff1 = 0
   
    if(!isTRUE(all.equal(aveEff1, 1))) return(0)
     

    PP =  matMulti1(PP, ginv(info.mat), C.ani, Z1.mat[ind[-length(ind)],]) 
            
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*%  PP %*% newX
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    #e.va = Re(svd(info.mat)$d)
    can.eff = e.va[-which(e.va < 1e-07)]/trt.rep
    len = length(can.eff)

    if(!isTRUE(all.equal(len/(nTrt-1), 1))) return(0)
    if(!isTRUE(all.equal((tr(PP) - len), newResDF))) return(0)  
   
    aveEff = 1/(mean(1/can.eff))
     
    return(aveEff)
}


obj.fun.RBD1 = function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {

  ans = try(try.obj.fun.RBD1(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, 
      blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF), silent = TRUE)
  
  ifelse(class(ans) == "try-error", 0, ans)
}



#####################################################################################
test.obj.fun.RBD = function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {
      #inforamtion in within runs and tags stratum while the cage information is eliminated
    info.mat = matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)
    PP = matMulti1(blk.proj, ginv(info.mat), C.cage,  Z1.mat[ind[-length(ind)],])
     
    if(any(abs(info.mat) > 1e-7)){     #check for any information in the cage stratum
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
    info.mat = matMulti(PP , Z1.mat[ind[-length(ind)],], C.ani)      
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    can.eff = e.va[-which(e.va < 1e-07)]/ani.Rep 

    aveEff1 = 1/mean(1/can.eff)
     
    if(is.na(aveEff1)) aveEff1 = 0
   
    if(!isTRUE(all.equal(aveEff1, 1))) return(0)
     

    PP =  matMulti1(PP, ginv(info.mat), C.ani, Z1.mat[ind[-length(ind)],]) 
           
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*%  PP %*% newX
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    #e.va = Re(svd(info.mat)$d)
    can.eff = e.va[-which(e.va < 1e-07)]/trt.rep
    len = length(can.eff)
     
    
    return(list(aveEff = 100/mean(1/can.eff), res = tr(PP) - len))
}


design.summary.RBD =
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
    print(summaryAovOnePhase(design.df, blk.str = "Cag/Ani", trt.str = "Trt"))

    cat("Phase 2 theoretical ANOVA:\n")
    print(summaryAovTwoPhase(design.df, blk.str2 = "Run", blk.str1 = "Cag/Ani", trt.str = "Tag + Trt"))

}
