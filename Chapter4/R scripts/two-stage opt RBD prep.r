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

newLETTERS = sort(sapply(strsplit(levels(interaction(LETTERS, LETTERS)), "\\."),
                function(x) paste(x, collapse = "", sep = "")))
                
#newLETTERS=LETTERS
sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R")
sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/optimTE/R")


is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


initialRBD = function(nTrt, bRep, nCag, tRep, nPlot) {
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

    print(summaryAovOnePhase(phase1DesignEX1, blk.str = "Cag/Ani", trt.str = "Trt"))


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
    print(phase1.mat)

    return(as.numeric(as.factor(t(phase1.mat[,,1]))))
}

#C++ matrix multiplication
code <-
  'arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat Z = Rcpp::as<arma::mat>(Z_);
 
  return wrap(Z*trans(Y) * X * Y * Z);'

matMulti <- cxxfunction( signature(X_="numeric",
Y_="numeric", Z_="numeric"),  
    ' arma::mat X = Rcpp::as<arma::mat>(X_);
      arma::mat Y = Rcpp::as<arma::mat>(Y_);
      arma::mat Z = Rcpp::as<arma::mat>(Z_);
 
      return wrap(Z*trans(Y) * X * Y * Z);', 
      plugin="RcppArmadillo")

#matMulti <- function(P, X, C) C %*% t(X) %*% P %*% X %*% C

code <-
  'arma::mat W = Rcpp::as<arma::mat>(W_);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat Z = Rcpp::as<arma::mat>(Z_);
 
  return wrap(W*Z*Y * X * Y * trans(Z)*W);'

matMulti1 <- cxxfunction( signature(W_="numeric",X_="numeric",
Y_="numeric", Z_="numeric"), 
'arma::mat W = Rcpp::as<arma::mat>(W_);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat Z = Rcpp::as<arma::mat>(Z_);
 
  return wrap(W*Z*Y * X * Y * trans(Z)*W);'
  , plugin="RcppArmadillo")

matMulti1 <- function(P, G, C, X)  P %*% X %*% C %*% G %*% C %*% t(X) %*% P 

################################################################################


initialBIBD = function(nTrt, bRep, nCag, tRep, nPlot) {
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt

    nAni = nTrt * bRep


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

print(summary.aov.onePhase(phase1DesignEX1, blk.str = "Cag/Ani", trt.str = "Trt"))

bibd = seq(1,nrow(phase1DesignEX1),nTrt+1)[1:nTrt]

temp = bibd

for(i in 2:ceiling(nCag/nTrt))
  temp = c(temp, bibd + temp[length(temp)])

if(sum(is.na(temp)) > 0)
  temp = temp[-which(is.na(temp))]

if(sum(temp>nAni) > 0)
  temp = temp[-which(temp>nAni)]


print(summary.aov.onePhase(phase1DesignEX1[-temp,], blk.str = "Cag/Ani", trt.str = "Trt"))

phase1DesignEX1 = phase1DesignEX1[-bibd,]

phase1DesignEX1$Ani =  newLETTERS[1:nrow(phase1DesignEX1)]

nAni = nrow(phase1DesignEX1)

n = nAni * tRep
k = nAni /nCag

aveEff.phase1 = (nTrt * (k-1))/(k * (nTrt-1))

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


#Parameter's of block structure of Phase 1
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
print(phase1.mat)

    return(as.numeric(as.factor(t(phase1.mat[,,1]))))
}



bestInd = bestInd2
#################################################################################

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
 
print(test(X.trt = mI(n), blk.proj = PP1, Rep = 1)[c("nCan", "can.eff")])


info.mat =  matMulti(blk.proj , Z1.mat[bestInd,], C.cage)

if(any(abs(info.mat) > 1e-7)){
PP = blk.proj - matMulti1(blk.proj, ginv( matMulti(blk.proj , Z1.mat[bestInd,], C.cage)), C.cage, Z1.mat[bestInd,]) 

PP1 = matMulti1(PP, ginv( matMulti(PP , Z1.mat[bestInd,], C.ani)), C.ani, Z1.mat[bestInd,]) 
} else{
PP1 = matMulti1( blk.proj, ginv( matMulti( blk.proj , Z1.mat[bestInd,], C.ani)), C.ani, Z1.mat[bestInd,]) 

} 

print(test(X.trt = X.trt[bestInd,], PP1, Rep=trt.rep))

print(fractions(test(X.trt = X.trt[bestInd,], PP1, Rep=trt.rep)$can))

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
            B = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                          as.character(phase1DesignEX1$Ani)),]$Cag,
            Ani = t(new.Z1.des),
            Trt = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                          as.character(phase1DesignEX1$Ani)),]$Trt)

print(summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt"))

levels(design.df$Ani) = rep(1:(nlevels(design.df$Ani)/nCag), nCag)

print(summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "B/Ani",
trt.str = "Tag + Trt"))


summaryAovOnePhase(design.df,  blk.str = "B/Ani", trt.str = "Trt", 
latex = TRUE, fixed.name = "\\tau")

summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "B/Ani",
trt.str = "Tag + Trt", latex = TRUE)

matrix(paste(as.numeric(design.df$B), LETTERS[design.df$Ani], sep = ""), 
        nrow = nBlk, ncol = nPlot, byrow = TRUE)


matrix(design.df$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE)

N = with(design.df, table(Trt,Run))
N %*% t(N)

N = with(design.df, table(Trt,Cag))
N %*% t(N)

(N = with(design.df, table(Cag,Run)))
N %*% t(N)

(N = with(design.df, table(Ani,Run)))
N %*% t(N)


save(design.df, file="RBDEX6.Rdata")

