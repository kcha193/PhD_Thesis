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

sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R")
sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/optimTE/R")

design.summary =
function(design.df, simple = TRUE) {


    n = nrow(design.df)
    nBlk = nlevels(design.df$Run)
    nPlot = nlevels(design.df$Tag)
    nCag = nlevels(design.df$B)

    nAni = nlevels(design.df$Ani)
    Run.mat = with(design.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
    Cag.mat = with(design.df, as.matrix(table(1:n, Cag)))
    Ani.mat = with(design.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
    Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))

    C.cage = (mI(nCag) - mK(nCag)) %x% mK(nAni)
    C.ani = mI(nCag) %x% (mI(nAni) - mK(nAni))
    cage.Rep = n/nCag

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
      print(test.RBD(X.trt = mI(n), blk.proj = PP1, Rep = 1)[c("nCan", "can.eff")])
  
      info.mat = matMulti( X = blk.proj, Y = Ani.mat, Z = C.cage)
  
      if (any(abs(info.mat) > 1e-07)) {
          PP = blk.proj - matMulti1(W = blk.proj, X =ginv(matMulti( X =blk.proj, Y = Ani.mat, Z = C.cage)),
              Y = C.cage, Z = Ani.mat)
  
          PP1 = matMulti1(W = PP, X = ginv(matMulti(X = PP, Y = Ani.mat, Z = C.ani)), Y = C.ani, Z = Ani.mat)
      } else {
          PP1 = matMulti1(W = blk.proj, X = ginv(matMulti(X = blk.proj, Y = Ani.mat, Z = C.ani)), Y = C.ani,
               Z = Ani.mat)
  
      }
  
  
      cat("Treatment design:\n")
      print(matrix(design.df$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE))
  
      cat("Treatment incidence matrix:\n")
      print((N = with(design.df, table(Trt, Run))))
  
      cat("Treatment concurrence matrix:\n")
      print(N %*% t(N))
  
      cat("Treatment efficiency:\n")
  
      print(test.RBD(X.trt = Trt.mat, PP1, Rep = n/ncol(Trt.mat)))
    }

    cat("Phase 1 theoretical ANOVA:\n")
    print(summaryAovOnePhase(design.df, blk.str = "Cag/Ani", trt.str = "Trt"))

    cat("Phase 2 theoretical ANOVA:\n")
    print(summaryAovTwoPhase(design.df, blk.str2 = "Run", blk.str1 = "Cag/Ani", trt.str = "Tag + Trt"))

}

design.df$Cag = as.factor(toupper(design.df$Cag))
#nCag   = nlevels(design.df$Cag)

#levels(design.df$Ani) = rep(1:(nlevels(design.df$Ani)/nCag), nCag)


design.summary(design.df)

summary.aov.onePhase(design.df, blk.str = "Cag/Ani", trt.str = "Tag + Trt")

summaryAovOnePhase(design.df, blk.str = "Run", trt.str = "B/Ani")

summary.aov.onePhase(design.df, blk.str = "Cag", trt.str = "Tag")




test(X.trt = Trt.mat, PP1, Rep = n/ncol(Trt.mat))$e.vec

trt.contr = test(X.trt = Trt.mat, PP1, Rep = n/ncol(Trt.mat))$e.vec

trt.contr1 = trt.contr[,1][as.numeric(design.df$Trt)]
trt.contr2 = trt.contr[,2][as.numeric(design.df$Trt)]


summaryAovTwoPhase(design.df, blk.str2 = "Run", blk.str1 = "Cag/Ani", 
trt.str = "Tag + Trt")



