sourceDir <- function(path, trace = TRUE, ...) {
  library(MASS)
  library(lattice)
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

newLETTERS = sort(sapply(strsplit(levels(interaction(LETTERS, LETTERS)), "\\."), function(x) paste(x, collapse = "", sep = "")))

sourceDir(path = "C:/Users/kcha193/My Dropbox/R functions for two-phase experiments")

LSD = function(x) {
    n = length(x)
    xx = (1:n) - 1
    d = numeric(n)

    for (i in 0:(n - 1)) d = rbind(d, x[((xx + i)%%n) + 1])

    return(d[-1, ])
}

   #phase1.mat[nBlk, , 1] = sapply(strsplit(diag(outer(rep(names(count),
            #    count), apply(phase1.mat[-nBlk, , 1], 2, function(x) paste(sort(x),
            #    collapse = ".", sep = "")), function(x, y) paste(x, y,
            #    sep = "."))), "\\."), function(x) x[1])

           # phase1.mat[nBlk, , 2] = as.character(phase1DesignEX1$Trt[match(phase1.mat[nBlk,
           #     , 1], phase1DesignEX1$Ani)])

nTrt = 5
bRep = 8
tRep = 2
nPlot = 8

test.design.even = function(nTrt, bRep, tRep, nPlot) {

    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt

    if ((nBlk%%1) > 0)
        stop("Total number of observation is not divisiable by the number of tag used!")

    nAni = nTrt * bRep

    c(nTrt, nAni, nPlot, n/nPlot)

    phase1DesignEX1 <- local({
        Ani = 1:nAni
        Trt = 1:nTrt

        data.frame(cbind(Ani, Trt))
    })
    phase1DesignEX1$Ani = as.factor(newLETTERS[phase1DesignEX1$Ani])
    phase1DesignEX1$Trt = as.factor(multiLetters(phase1DesignEX1$Trt))

    print(summary.aov.onePhase(phase1DesignEX1, blk.str = "Ani", trt.str = "Trt"))

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
            trt = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani ==
                names(count)[i])])
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

            phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, ani.limit,
                trt.limit)

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
           for(j in 1:nPlot){
            phase1.mat[nBlk, j, 1] = fillIn(nBlk, j, phase1.mat, count,
                    rep(Z1.rep, 2), c(nBlk/nTrt, trt.rep))
           
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

                  phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count,
                    rep(Z1.rep, 2), c(nBlk/nTrt, nTrt))

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
                  phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count,
                    ani.limit, c(nBlk/nTrt, trt.rep/tRep))

                  phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani ==
                    phase1.mat[i, j, 1])])

                  count[phase1.mat[i, j, 1]] = count[phase1.mat[i, j, 1]] -
                    1

                  # count = sort(count, d=TRUE)
                  if (any(count == 0))
                    count = count[-which(count == 0)]
                }
            }

        }
    }
    print(phase1.mat)

    Z1.des = as.factor(t(phase1.mat[, , 1]))
    trt.des = as.factor(t(phase1.mat[, , 2]))

    design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

    design.df = cbind(design.df, Ani = Z1.des, Trt = trt.des)


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
    cat("Trt:", ncol(Trt.mat), "Ani:", ncol(Ani.mat), "Tag:", ncol(Tag.mat),
        "Run:", ncol(Run.mat), "\n")

    cat("Animal efficiency:\n")
    print(test(X.trt = Ani.mat, blk.proj, Rep = n/ncol(Ani.mat)))
    
    cat("Treatment efficiency:\n")
    print(test(X.trt = Trt.mat, (mI(n) - Pb) %*% (mI(n) - Pb1), Rep = n/ncol(Trt.mat)))
    
    summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani",
        trt.str = "Tag + Trt")

}

test.design.even(2, 2, 2, 4)
test.design.even(2, 4, 2, 4)
test.design.even(2, 6, 2, 4)
test.design.even(2, 8, 2, 4)

test.design.even(2, 3, 2, 4)
test.design.even(2, 5, 2, 4)
test.design.even(2, 7, 2, 4)
test.design.even(2, 9, 2, 4)
test.design.even(2, 11, 2, 4)

test.design.even(2, 2, 4, 4)
test.design.even(2, 3, 4, 4)  #<-
test.design.even(2, 4, 4, 4)
test.design.even(2, 5, 4, 4) #<-
test.design.even(2, 6, 4, 4)
test.design.even(2, 7, 4, 4)  #<-
test.design.even(2, 8, 4, 4)
test.design.even(2, 9, 4, 4)  #<-
test.design.even(2, 10, 4, 4)  
test.design.even(2, 11, 4, 4)  

test.design.even(2, 4, 2, 8)
test.design.even(2, 6, 2, 8)
test.design.even(2, 8, 2, 8)
test.design.even(2, 10, 2, 8)
test.design.even(2, 12, 2, 8)
test.design.even(2, 14, 2, 8)

test.design.even(2, 2, 4, 8)
test.design.even(2, 3, 4, 8)

test.design.even(2, 4, 4, 8)
test.design.even(2, 5, 4, 8)
test.design.even(2, 6, 4, 8)
test.design.even(2, 7, 4, 8)
test.design.even(2, 8, 4, 8)
test.design.even(2, 9, 4, 8)
test.design.even(2, 10, 4, 8)
test.design.even(2, 11, 4, 8)

test.design.even(4, 3, 2, 4)
test.design.even(4, 4, 2, 4)
test.design.even(4, 5, 2, 8)
test.design.even(4, 6, 2, 4)
test.design.even(4, 7, 2, 4)

test.design.even(6, 2, 2, 4)
test.design.even(6, 4, 2, 4)
test.design.even(6, 6, 2, 4)
test.design.even(6, 8, 2, 4)
test.design.even(6, 10, 2, 4)

test.design.even(6, 3, 2, 4)

test.design.even(8, 2, 2, 4)

test.design.even(10, 2, 2, 4)


test.design.even(3, 2, 2, 4)
test.design.even(3, 4, 2, 4)
test.design.even(3, 6, 2, 4)
test.design.even(3, 8, 2, 4)
test.design.even(3, 10, 2, 4)

test.design.even(5, 8, 2, 8)

test.design.even(6, 4, 2, 8)

test.design.even(7, 2, 2, 4)
test.design.even(7, 4, 2, 4)
test.design.even(7, 8, 2, 4)
test.design.even(7, 10, 2, 4)

#nTrt, bRep, tRep, nPlot
test.design.even(7, 4, 2, 8)
test.design.even(7, 4, 4, 8)

test.design.odd = function(nTrt, bRep, tRep, nPlot) {
    
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    
    nAni = nTrt * bRep
    
    
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
    
    fillIn = function(run, tag, phase1.mat, count, ani.limit, trt.limit) {
        
        for (i in 1:length(count)) {
            
            check1 = sum(phase1.mat[run, , 1] == names(count)[i]) < ani.limit[2]
            check2 = sum(phase1.mat[, tag, 1] == names(count)[i]) < ani.limit[1]
            trt = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == 
                names(count)[i])])
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
    
    if (nZ1 < 8) {
        x = sort(seq(1, nZ1, 8))
    } else {
        x = sort(c(seq(1, nZ1, 8), seq(6, nZ1, 8)))
    }
    
    count.ind = t(sapply(x, function(y) y:(y + 2)))
    count.total = matrix(1:nZ1, nrow = nBlk/tRep, ncol = nPlot, byrow = TRUE)
    
    ani.limit = c(nBlk, nPlot)/nZ1
    trt.limit = c(nBlk, nPlot)/nTrt
    
    phase1.mat = array("-1", c(nBlk, nPlot, 2))
    
    len = nBlk
    
    for (k in 1:(nBlk/tRep)) {
        count1 = count[count.ind[k, ]]
        countT = count.total[k, ]
        
        for (l in 1:(nBlk/(nBlk/tRep))) {
            i = l + (k - 1) * (nBlk/(nBlk/tRep))
            
            for (j in 1:(nPlot - 1)) {
                
                phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count1, ani.limit, 
                  trt.limit)
                
                phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == 
                  phase1.mat[i, j, 1])])
                
                count1[phase1.mat[i, j, 1]] = count1[phase1.mat[i, j, 1]] - 
                  1
                
                count1 = sort(count1, d = TRUE)
                if (any(count1 == 0)) 
                  count1 = count1[-which(count1 == 0)]
            }
            
            last = countT[-match(count.ind[k, ], countT)]
            
            phase1.mat[i, nPlot, 1] = names(count)[last]
            phase1.mat[i, nPlot, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == 
                phase1.mat[i, nPlot, 1])])
            # count[phase1.mat[i, nPlot, 1]] = count[phase1.mat[i, nPlot, 1]] - 1
        }
    }
    
    
    print(phase1.mat)
    
    Z1.des = as.factor(t(phase1.mat[, , 1]))
    trt.des = as.factor(t(phase1.mat[, , 2]))
    
    design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))
    
    design.df = cbind(design.df, Ani = Z1.des, Trt = trt.des)
    
    
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
    cat("Trt:", ncol(Trt.mat), "Ani:", ncol(Ani.mat), "Tag:", ncol(Tag.mat), 
        "Run:", ncol(Run.mat), "\n")
    
    cat("Animal efficiency:\n")
    print(test(X.trt = Ani.mat, blk.proj, Rep = n/ncol(Ani.mat)))
    
    summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")
    
}


# nTrt, bRep, tRep, nPlot
test.design.odd(2, 2, 3, 4)
test.design.odd(2, 4, 3, 4)
test.design.odd(2, 6, 3, 4)
test.design.odd(2, 8, 3, 4)
test.design.odd(2, 10, 3, 4)
test.design.odd(2, 12, 3, 4)


test.design.odd(2, 4, 3, 8)
test.design.odd(2, 8, 3, 8)
test.design.odd(2, 12, 3, 8)


test.design.even(2, 8, 3, 6)


test.design.odd(4, 2, 3, 4)





    
 