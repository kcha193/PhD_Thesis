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

sourceDir(path = "C:/Users/kcha193/My Dropbox/R functions for two-phase experiments")


design.summary = function(design.df, trtFirst = TRUE) {
    
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
    
    
    cat("Animal design:\n")
    print(matrix(design.df$Ani, nrow = nBlk, ncol = nPlot, byrow = TRUE))
    
    cat("Animal incidence matrix:\n")
    print((N = with(design.df, table(Ani, Run))))
    
    cat("Animal concurrence matrix:\n")
    print(N %*% t(N))
    
    cat("Animal efficiency:\n")
    print(test(X.trt = Ani.mat, blk.proj, Rep = n/ncol(Ani.mat)))
    
    cat("Treatment design:\n")
    print(matrix(design.df$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE))
    
    cat("Treatment incidence matrix:\n")
    print((N = with(design.df, table(Trt, Run))))
    
    cat("Treatment concurrence matrix:\n")
    print(N %*% t(N))
    
    cat("Treatment efficiency:\n")
    print(test(X.trt = Trt.mat, (mI(n) - Pb) %*% (mI(n) - Pb1), Rep = n/ncol(Trt.mat)))
    
    # blk.proj = (mI(n) - Pb)
    #PP = blk.proj %*% Ani.mat %*% (mI(ncol(Ani.mat)) - mK(ncol(Ani.mat))) %*% invInfMat(blk.proj, 
    #    Ani.mat, mI(ncol(Ani.mat)) - mK(ncol(Ani.mat))) %*% (mI(ncol(Ani.mat)) - mK(ncol(Ani.mat))) %*% 
    #     t(Ani.mat) %*% blk.proj
    
    #cat("Treatment efficiency:\n")
    #print(test(X.trt = Trt.mat, blk.proj = PP, Rep = n/ncol(Trt.mat)))
    
    cat("Phase 1 theoretical ANOVA:\n")
    print(summary.aov.onePhase(design.df, blk.str = "Ani", trt.str = "Trt"))

    cat("Phase 2 theoretical ANOVA:\n")
    summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")

    
} 

design.summary(design.df)


design.summary(design.df, FALSE)



summary.aov.onePhase(design.df, blk.str = "Run", trt.str = "Ani")
   Cag = LETTERS[1:3]
   Ani = 1:6
    
  levels(design.df$Ani) = sort(levels(interaction(Cag, Ani)))

 design.df$Cag = factor(sapply(strsplit(as.character(design.df$Ani), "\\."), function(x) x[1]))
  design.df$Ani = factor(sapply(strsplit(as.character(design.df$Ani), "\\."), function(x) x[2]))




  trt.contr = test(X.trt = Trt.mat, (mI(n) - Pb) %*% (mI(n) - Pb1), Rep = n/ncol(Trt.mat))$e.vec

trt.contr1 = trt.contr[,1][as.numeric(design.df$Trt)]
trt.contr2 = trt.contr[,2][as.numeric(design.df$Trt)]
trt.contr3 = trt.contr[,3][as.numeric(design.df$Trt)]
trt.contr4 = trt.contr[,4][as.numeric(design.df$Trt)]
trt.contr5 = trt.contr[,5][as.numeric(design.df$Trt)]


summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt",trt.contr = list(Tag = NA,
                                       Trt = list(  "1" = trt.contr1,
                                                    "2" = trt.contr2,
                                                    "3" = trt.contr3,
                                                    "4" = trt.contr4,
                                                    "5" = trt.contr5)))
summary.aov.onePhase(design.df, blk.str = "Ani", trt.str = "Trt", latex = TRUE)
  
  
  summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt", latex = TRUE)  
  
   summary.aov.twoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag*Trt") 
   
   