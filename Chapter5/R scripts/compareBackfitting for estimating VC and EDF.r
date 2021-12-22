

resProj = function(design1.df, true.VC) {
    
    n = nrow(design1.df)    
    nBlk = nlevels(design1.df$Run)
    nPlot = nlevels(design1.df$Tag)
       
    nAni = nlevels(design1.df$Ani)
    Run.mat = with(design1.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design1.df, as.matrix(table(1:n, Tag)))
    Ani.mat = with(design1.df, as.matrix(table(1:n, Ani)))
    Trt.mat = with(design1.df, as.matrix(table(1:n, Trt)))
    
    browser()
    
    R1 = Run.mat %*% t(Run.mat)
    D1 = R1 * true.VC[3]
        
    R2 = Ani.mat %*% t(Ani.mat)
    D2 = R2 * true.VC[2]
        
    V = diag(1, n) * true.VC[1] + D1 +  D2
  
    C.ani =  (identityMat(nAni) - K(nAni))
    
    Pb = projMat(Run.mat)
    Pb1 = projMat(Tag.mat)
      
    blk.proj = (identityMat(n) - Pb) %*% (identityMat(n) - Pb1)
    
    PP =  matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani, Ani.mat)
    
    X = Trt.mat    
    
    resProj = PP - PP %*% X %*%  ginv(t(X) %*% PP %*% X ) %*% t(X) %*% PP
     
    return( ((tr(resProj %*% V))^2)/tr(resProj %*% V %*% resProj %*% V) )
}


design2 = design
















summaryAovTwoPhase(design, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")

REML.EDF = numeric(nS)
LC.EDF = numeric(nS)
REAL.EDF = numeric(nS)

REML.VC = numeric(nVc)
LC.VC = numeric(nVc)
REAL.VC = numeric(nVc)

gamma.run = 10
gamma.ani = 5
VC.resid = 1
nSim = 1000

pb <- txtProgressBar(min = 0, max = nSim, style = 3)
iSim = 0
# simulation
while (iSim < nSim) {
    setTxtProgressBar(pb, iSim + 1)
    
    # Introduce the variation
    run.eff = rnorm(nlevels(design$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
    ani.eff = rnorm(nlevels(design$Ani), mean = 0, sd = sqrt(gamma.ani * VC.resid))
    trt.eff = c(1,3) 
    tag.eff = runif(nlevels(design$Tag), 0, 1)
    res.eff = rnorm(nrow(design), mean = 0, sd = sqrt(VC.resid))
    gm = 10
    
    (real.VC = c(VC.resid, (gamma.ani * VC.resid), (gamma.run * VC.resid)))
    
    y = gm + with(design, run.eff[Run] + ani.eff[Ani] + tag.eff[Tag] + trt.eff[Trt]) + res.eff
    
    
    # Using rubin's method to
    emY = y
    
    # Compute the parameters using the theoritical ANOVA table
 
     n = length(y)
    design.df = design
    
    nAni = nlevels(design.df$Ani)
    Run.mat = with(design.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
#    Cag.mat = with(design.df, as.matrix(table(1:n, Cag)))
    Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
    Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))

    # Treatment design matrix
    X = cbind(1, Trt.mat, Tag.mat)
    
    L = X %*% ginv(t(X) %*% X) %*% t(X)
    
    # Run design matrix
    Z1 = Run.mat
      
    # Animal design matrix
    Z2 = Ani.mat
     
    # Least square estimates
    thetaHat = ginv(t(X) %*% X) %*% t(X) %*% emY
    
    #thetaHat =  trtProj(design.df) %*% emY
     
    yc = emY - X %*% thetaHat
    
    b2Hat = ginv(t(Z2) %*% Z2 ) %*% t(Z2) %*% yc
       
    conv = matrix(0, nrow = 2, ncol = 3)
    
    conv[nrow(conv), ] = var.comp = real.VC
    
    sigma.run = var.comp[3]
    sigma.ani = var.comp[2]
    sigma.base = var.comp[1]
   
   
    counter = 0
    # print(var.comp) Check over 100 iterations
    while (!all(apply(conv, 1, function(x) isTRUE(all.equal(x, conv[nrow(conv), ], tolerance = 1e-10))))) {
        
        sigma.run = var.comp[3]
        sigma.ani = var.comp[2]
        sigma.base = var.comp[1]
                
        R1 = diag(1, ncol(Run.mat))
        D1 = R1 * sigma.run
        
         R2 = diag(1, ncol(Ani.mat))
        D2 = R2 * sigma.ani
        
        A = diag(1, n)
        # A[miss,]<-0
        sig = A * sigma.base
        
        yc = emY - X %*% thetaHat - Z2 %*% b2Hat
  
        # Estimate run effects
        b1Hat = ginv(t(Z1) %*% ginv(sig) %*% Z1 + ginv(D1)) %*% t(Z1) %*% ginv(sig) %*% yc
        
        yc = emY - X %*% thetaHat  - Z1 %*% b1Hat
        
        # Estimate animal effects without run effects and fixed effects 
        b2Hat = ginv(t(Z2) %*% ginv(sig) %*% Z2 + ginv(D2)) %*% t(Z2) %*% ginv(sig) %*% yc
        
        yc = emY - Z2 %*% b2Hat - Z1 %*% b1Hat
        
        # Estimate theta
        thetaHat = ginv(t(X) %*% ginv(sig) %*% X) %*% t(X) %*% ginv(sig) %*% yc
        
        #yc = emY - X %*% thetaHat - Z2 %*% b2Hat
   
        # Estimate the variance componenets
        e = resid(lm(emY ~ Run + Ani + Trt + Tag, design))
        
        #sigma.base = as.numeric((1/(11)) * (t(e) %*% A %*% e + 
        #        sum(diag(ginv(t(Z1) %*% A %*% Z1/sigma.base + R1/sigma.run) %*% t(Z1) %*% A %*% Z1)) + 
        #        sum(diag(ginv(t(Z2) %*% A %*% Z2/sigma.base + R2/sigma.ani) %*% t(Z2) %*% A %*% Z2))))

        
        sigma.base =  as.numeric((t(e) %*% A %*% e ) /(11 -  
         sum(diag(ginv(t(Z1) %*% A %*% Z1 + R1* sigma.base/sigma.run) %*% t(Z1) %*% A %*% Z1))-
         sum(diag(ginv(t(Z2) %*% A %*% Z2 + R2* sigma.base/sigma.ani) %*% t(Z2) %*% A %*% Z2))))
        
       
        sigma.run = as.numeric((1/
                sum(diag(ginv(t(Z1) %*% A %*% Z1 + R1* sigma.base/sigma.run) %*% t(Z1) %*% A %*% Z1))) * 
                (t(b1Hat) %*% R1 %*% b1Hat + 
              sum(diag(ginv((t(Z1) %*% A %*% Z1)/sigma.base + R1/sigma.run) %*% R1))))
        
        sigma.ani = as.numeric((1/
                sum(diag(ginv(t(Z2) %*% A %*% Z2 + R2* sigma.base/sigma.ani) %*% t(Z2) %*% A %*% Z2))) * 
                (t(b2Hat) %*% R2 %*% b2Hat + 
              sum(diag(ginv((t(Z2) %*% A %*% Z2)/sigma.base + R2/sigma.ani) %*% R2))))
        
            
        var.comp = c(sigma.base, sigma.ani, sigma.run)
        
        print(var.comp)
        
        conv = rbind(conv, t(var.comp))
        conv = conv[-1, ]
        counter = counter + 1
        if (counter > 1000) {
            break
        }
        
    }
    
     aov.table = summaryAovTwoPhase(design, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt", response = y)
    
    
    tmp = try(getVcEDF(aov.table = aov.table, row.MS = NA, true.VC = real.VC, neg.VC = FALSE), TRUE)
          #try(getVcEDF(aov.table = aov.table, row.MS = NA, true.VC = var.comp, neg.VC = TRUE), TRUE)
    
    
    # tmp
    
    if (class(tmp) == "try-error") {
        iSim = iSim - 1
        next
    }
   
       tmpS = tmp$S
        tmpV = tmp$V

    #REML.EDF = rbind(REML.EDF, tmpS[, "REML.EDF"])
    #LC.EDF = rbind(LC.EDF, tmpS[, "LC.EDF"])
    #REAL.EDF = rbind(REAL.EDF, tmpS[, "TRUE.EDF"])
    
    REML.VC = rbind(REML.VC, tmpV[, "REML.var.comp"])
    LC.VC = rbind(LC.VC, tmpV[, "LC.var.comp"])
    REAL.VC = rbind(REAL.VC, var.comp)
    
    iSim = iSim + 1
}
close(pb) 



getVcEDF(aov.table = aov.table, row.MS = NA, true.VC = real.VC, neg.VC = neg.VC)


     apply(REML.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
     apply(LC.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
      apply(REAL.EDF[-1,], 2, function(x) median(x, na.rm = TRUE))
      apply(REAL.EDF[-1,], 2, function(x) mean(x, na.rm = TRUE))
    
     apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
    apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))
     apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE))


     plot(REML.VC[-1,1], LC.VC[-1,1])
       plot(REML.VC[-1,2], LC.VC[-1,2])
       plot(REML.VC[-1,3], LC.VC[-1,3])
     
        plot(REML.VC[-1,1], REAL.VC[-1,1])
         plot(REML.VC[-1,2], REAL.VC[-1,2])
          plot(REML.VC[-1,3], REAL.VC[-1,3])
   
     
     
     