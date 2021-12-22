initial = initialRBD.nolimit

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

Z1.mat = matrix(0, ncol = nAni, nrow = nAni * Z1.rep)
Z1.mat[cbind(1:(nAni * Z1.rep), Z1.des)] = 1
Z1.rep = nrow(Z1.mat)/ncol(Z1.mat)

new.Z1.mat=  Z1.mat

trt.rep = n/nTrt

trt.des = as.numeric(phase1DesignEX1$Trt[match(newLETTERS[Z1.des], phase1DesignEX1$Ani)])
X.trt = matrix(0, ncol = nTrt, nrow = nTrt * trt.rep)
X.trt[cbind(1:(nTrt * trt.rep), trt.des)] = 1

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

matrix(as.numeric(design.df$Cag), nrow = nBlk, ncol = nPlot, byrow = TRUE)
matrix(design.df$Ani, nrow = nBlk, ncol = nPlot, byrow = TRUE)

matrix(design.df$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE)




threeStageRBD.spe = function(init, iter = 50000, obj.fun = obj.fun.new) {
    print(c(nTrt, bRep, nCag, tRep, nPlot))
   
    newInit = c(init, init[1])

    old = obj.fun(newInit)

    cat("Level: 1, Finding the temperatures for both stages\n")

    
    temp = initialTemp(iter, newInit, swap.stage2.new, obj.fun)

    res.stage2 = temp$newInit
    y2.stage2 = temp$y2
    y2low.stage2 = temp$y2low
    
    cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
 
    newInit = res.stage2
    
    cur = obj.fun(newInit)
    print(test.obj.fun.new(newInit))
    
    while ((cur - old) > 1e-06) {
        
        temp = initialTemp(iter, newInit, swap.stage2.new, obj.fun)
    
        res.stage2 = temp$newInit
        y2.stage2 = temp$y2
        y2low.stage2 = temp$y2low
                
        newInit = res.stage2

        old = cur
        
        cur = obj.fun(newInit)
        print(test.obj.fun.new(newInit))     
   
        cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
    
     }
      

    old = obj.fun(newInit)
  
    cat("Level: 2, Current temp: ", y2.stage2, "\n")

    res <- optim( newInit, obj.fun, swap.stage2.new, method = "SANN", 
                  control = list(maxit = iter, temp = y2.stage2,
                  tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

    cur = obj.fun(res$par)

    print(test.obj.fun.new(res$par))
    
    if(cur == 1) return(res$par[-length(res$par)])

    while ((cur - old) > 1e-06) {
        res <- optim(res$par, obj.fun, swap.stage2.new, method = "SANN", 
                      control = list(maxit = iter, temp = y2.stage2,
                      tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))

        old = cur
        cur = obj.fun(res$par)

        print(test.obj.fun.new(res$par))


        if(cur == 1) return(res$par[-length(res$par)])
    }



    inter.stage2 = exp(log(y2.stage2/y2low.stage2)/8)

    i = 1
    unstable = FALSE

    # for (i in 1:8) {
    while (i < 9 || unstable) {
        unstable = FALSE

        y2.stage2 = y2.stage2/inter.stage2

        
        old = obj.fun(res$par)

        cat("Level: ", i + 2, ", Current temp: ", y2.stage2, "\n")
        res <- optim(res$par, obj.fun, swap.stage2.new, method = "SANN", 
                      control = list(maxit = iter, temp = y2.stage2,
                      tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))


        cur = obj.fun(res$par)
        print(test.obj.fun.new(res$par))

        if(cur == 1) return(res$par[-length(res$par)])

        while ((cur - old) > 1e-06) {
            unstable = TRUE
            res <- optim(res$par, obj.fun, swap.stage2.new, method = "SANN", 
                          control = list(maxit = iter, temp = y2.stage2, 
                          tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1))
    
            old = cur
            cur = obj.fun(res$par)
            print(test.obj.fun.new(res$par))

            if(cur == 1) return(res$par[-length(res$par)])      
        }
                                                                  
        i = i + 1
    }

    print(c(nTrt, bRep, nCag, tRep, nPlot))

    return(res$par[-length(res$par)])
}




design.df$Ani =  as.factor(c(1,2,3,4,5,6,7,1,
                   2,1,4,3,6,5,1,7,
                   2,3,4,5,6,7,1,2,
                   3,2,5,4,7,6,2,1,
                   3,4,5,6,7,1,2,3,
                   4,3,6,5,1,7,3,2,
                   4,5,6,7,1,2,3,4,
                   5,4,7,6,2,1,4,3,
                   5,6,7,1,2,3,4,5,
                   6,5,1,7,3,2,5,4,
                   6,7,1,2,3,4,5,6,
                   7,6,2,1,4,3,6,5,
                   7,1,2,3,4,5,6,7,
                   1,7,3,2,5,4,7,6))
                   
design.df$Trt =  as.factor(letters[c(1,2,3,4,5,6,7,1,
                   2,1,4,3,6,5,1,7,
                   2,3,4,5,6,7,1,2,
                   3,2,5,4,7,6,2,1,
                   3,4,5,6,7,1,2,3,
                   4,3,6,5,1,7,3,2,
                   4,5,6,7,1,2,3,4,
                   5,4,7,6,2,1,4,3,
                   5,6,7,1,2,3,4,5,
                   6,5,1,7,3,2,5,4,
                   6,7,1,2,3,4,5,6,
                   7,6,2,1,4,3,6,5,
                   7,1,2,3,4,5,6,7,
                   1,7,3,2,5,4,7,6)])
          
