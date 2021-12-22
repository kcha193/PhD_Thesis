#Compare the Cage allocations

nTrt = 3; 
bRep = 4; 
nCag = 4; 
tRep = 2; 
nPlot = 4
 
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

betRun = Pb - K(n)
betTag = Pb1 - K(n)
withBlock = (identityMat(n) - Pb) %*% (identityMat(n) - Pb1)

blk.proj = withBlock

# Parameter's of block structure of Phase 1
nZ1 = nrow(phase1DesignEX1)
Z1.rep = n/nZ1

Z1.des.sortRun = as.numeric(initialRBD(nTrt = nTrt, bRep = bRep, nCag = nCag, 
                    tRep = tRep, nPlot = nPlot))

Z1.des.sortRun = phase1DesignEX1[match(multiLetters(Z1.des.sortRun), as.character(phase1DesignEX1$Ani)), ]

levels(Z1.des.sortRun$Ani) = multiLetters(rep(1:(nlevels(Z1.des.sortRun$Ani)/nCag), nCag) )

matrix(paste(Z1.des.sortRun$Cag, Z1.des.sortRun$Ani, sep = ""), nrow = nBlk, ncol = nPlot, byrow = TRUE) 
matrix(Z1.des.sortRun$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE) 

Z1.des.sortTag = as.numeric(initialRBD.nolimit1(nTrt = nTrt, bRep = bRep, nCag = nCag, 
            tRep = tRep, nPlot = nBlk))
            
Z1.des.sortTag = phase1DesignEX1[match(multiLetters(Z1.des.sortTag), as.character(phase1DesignEX1$Ani)),]

levels(Z1.des.sortTag$Ani) = multiLetters(rep(1:(nlevels(Z1.des.sortTag$Ani)/nCag), nCag) )

matrix(paste(Z1.des.sortTag$Cag, Z1.des.sortTag$Ani, sep = ""), nrow = nBlk, ncol = nPlot, byrow = TRUE) 
matrix(Z1.des.sortTag$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE) 

################################################################################
 
design.df = optCRD(nTrt  = 4, bRep  = 2, tRep  = 2, nPlot = 4, iter  = 1000)

design.summary.CRD(design.df, FALSE)
 
n = nrow(design.df)
nBlk = nlevels(design.df$Run)
nPlot = nlevels(design.df$Tag)
Run.mat = with(design.df, as.matrix(table(1:n, Run)))
Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))
Pb = projMat(Run.mat)
Pb1 = projMat(Tag.mat)

(ave.eff = test.CRD(X.trt = Trt.mat, (diag(n) - Pb) %*% (diag(n) - Pb1), Rep = n/ncol(Trt.mat))$ave.eff )

design1.df = optRBD( nTrt = 4, bRep  = 2, nCag  = 2, tRep  = 2, nPlot = 4, iter  = 3000, 
upperValue = ave.eff, confoundCag = FALSE)

design.summary.RBD(design1.df, FALSE)
design.summary.RBD(design1.df, FALSE, latex = TRUE)
 
design2.df = optRBD( nTrt = 4, bRep  = 2, nCag  = 2, tRep  = 2, nPlot = 4, iter  = 3000, 
upperValue = ave.eff, confoundCag = TRUE)
 
design.summary.RBD(design2.df, FALSE)
 
 
 
 
 
 
 
 
 
 