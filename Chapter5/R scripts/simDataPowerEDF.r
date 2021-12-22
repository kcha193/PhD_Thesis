
library(snowfall)  #Use parallel programming to speed up the computation by 4 times
sfInit(parallel = TRUE, cpus = 4)

sourceDir <- function(path, trace = TRUE, ...) {
    library(MASS)
    
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        if (trace) 
            cat(nm, ":")
        source(file.path(path, nm), ...)
        if (trace) 
            cat("\n")
    }
}

design.df = optCRD(nTrt = 2, bRep  = 3, tRep  = 2, nPlot = 4, iter  = 1000)

(anovaTable = summaryAovTwoPhase(design.df = design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")$A)



simDataPowerEDF = function(design.df, trt.Eff, gamma.R, gamma.A, VC.resid, nSim, row.MS = NA, neg.VC = TRUE) {
    
    sfInit(parallel = TRUE, cpus = 4)
    
    temp = data.frame()
    
    for (i in 1:nrow(trt.Eff)) {
        cat("trt.Eff ", i, "\n")
        
        trt.eff = trt.Eff[i, ]
        
        for (k in 1:length(gamma.R)) {
            gamma.run = gamma.R[k]
            cat("gamma.run ", gamma.run, "\n")
            
            
            test = function(j) {
                # sourceDir(path = '/home/kcha193/infoDecompuTE/R', trace = FALSE)
                sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R", trace = FALSE)
                
                o.MS = t(as.matrix(rep(0, 3)))
                r.MS = t(as.matrix(rep(0, 3)))
                
                counter = 0
                # simulation
                
                while (counter < nSim) {
                  
                  run.eff = rnorm(nlevels(design.df$Run), mean = 0, sd = sqrt(gamma.run * VC.resid))
                  ani.eff = rnorm(nlevels(design.df$Ani), mean = 0, sd = sqrt(gamma.A[j] * VC.resid))
                  tag.eff = runif(nlevels(design.df$Tag), 0, 1)
                  res.eff = rnorm(nrow(design.df), mean = 0, sd = sqrt(VC.resid))
                  gm = 10
                  
                  y = gm + with(design.df, run.eff[Run] + ani.eff[Ani] + tag.eff[Tag] + trt.eff[Trt]) + res.eff
                  
                  
                  aov.table = summaryAovTwoPhase(design.df = design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt", 
                    response = y)
                  
                  tmp = try(getVcEDF(aov.table = aov.table, true.VC = c(VC.resid, (gamma.A[j] * VC.resid), (gamma.run * VC.resid)), 
                    neg.VC = neg.VC), TRUE)
                  
                  if (class(tmp) == "try-error") {
                    counter = counter + 1
                    next
                  }
                  
                  o.MS = rbind(o.MS, c(2, as.numeric(aov.table$A[7:8, "MS"])))
                  
                  EDF = max(tmp$S["Within RunBetweenAniResidual", c("REML.EDF", "LC.EDF")], na.rm = TRUE)
                  
                  r.VC = tmp$V[, which(tmp$S["Within RunBetweenAniResidual", c("REML.EDF", "LC.EDF")] == EDF)]
                  
                  r.MS = rbind(r.MS, c(EDF, as.numeric(aov.table$A[7, "MS"]), sum(as.numeric(aov.table$A[8, 2:4]) * r.VC)))
                  
                  counter = counter + 1
                }
                
                o.MS = o.MS[-1, ]
                r.MS = r.MS[-1, ]
                
                return(c(paste(trt.eff, sep = "", collapse = ","), c(gamma.run, gamma.A[j], mean(o.MS[, 1]), 
                      sum((1 - pf(o.MS[, 2]/o.MS[, 3], 1, o.MS[, 1])) < 0.05)/nSim, mean(o.MS[, 3]),
                         mean(r.MS[, 1]), 
                         sum((1 - pf(r.MS[, 2]/r.MS[, 3], 1, r.MS[, 1])) < 0.05)/nSim, mean(r.MS[, 3]))))
            }
            
            sfExport("gamma.A", "gamma.run", "trt.eff", "VC.resid", "design.df", "sourceDir", "nSim", "getVcEDF", "fracToNum", 
                "extractName")
            
            sfLibrary(MASS)
            
            # old = proc.time()
            temp = rbind(temp, t(sfSapply(1:length(gamma.A), test)))
            # proc.time() - old print(dim(temp))
        }
        
    }
    sfStop()
    
    colnames(temp) = c("trt.eff", "gamma.run", "gamma.ani", "DF", "Initial Power", "MS", "EDF", "Adjusted Power", "Adjusted MS")
    
    temp$DF = as.numeric(as.character(temp$DF))
    temp$"Initial Power" = as.numeric(as.character(temp$"Initial Power"))
    temp$EDF = as.numeric(as.character(temp$EDF))
    temp$"Adjusted Power" = as.numeric(as.character(temp$"Adjusted Power"))
    temp$gamma.ani = as.numeric(as.character(temp$gamma.ani))
    temp
} 


trt.Eff = t( as.matrix( c(0, 4)))

old = proc.time()
adjustZero = 
    simDataPowerEDF(design.df = design.df, 
                    trt.Eff = trt.Eff,
                    gamma.R =  c(0, 0.25, 1, 4, 100), 
                    gamma.A = c(10^((-8:8)/2)), 
                    VC.resid = 1, 
                    nSim = 10000, 
                    neg.VC = TRUE)
    sfStop()
 proc.time() - old



old = proc.time()
remainNeg = 
    simDataPowerEDF(design.df = design.df, 
                    trt.Eff = trt.Eff,
                    gamma.R =  c(0, 0.25, 1, 4, 100), 
                    gamma.A = c(10^((-8:8)/2)), 
                    VC.resid = 1, 
                    nSim = 10000, 
                    neg.VC = FALSE)
    sfStop()
 proc.time() - old



old = proc.time()

sim4Tags = simDataChisqP(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100), 
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1057, neg.VC = FALSE, checkFun = FALSE)
 proc.time() - old

edfPlot(sim4Tags)

(aov.table = summaryAovTwoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani",
trt.str = "Tag + Trt", response = rnorm(nrow(design.df)))$A)

(rowNames = extractName(rownames(aov.table)))

xx =  simDataChisqP(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10000, seed = 1, 
            neg.VC = FALSE,
            row.MS = rowNames[c(8, 11)])

sim4Tags$tempS = rbind(sim4Tags$tempS, xx$tempS)
 
edfPlot(sim4Tags)

sim4Tagstemp = sim4Tags
 
sim4Tagstemp$tempS = sim4Tagstemp$tempS[,-5]
 
pdf("CRD232Tag4Unadjusted.pdf", width = 12)
edfPlot(sim4Tagstemp)
dev.off()

old = proc.time()
sim4Tags1 = simDataChisqP(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = c(0, 0.25, 1, 4, 100),
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 4, nVc = 3, nSim = 10000, seed = 1057, neg.VC = TRUE, checkFun = FALSE)
proc.time() - old
            
edfPlot(sim4Tags1)


xx =  simDataChisqP(design =  design.df,
            blk.str2 = "Run", blk.str1 = "Ani",
            trt.str = "Tag + Trt",
            gamma.R = NA,
            gamma.A =c(10^((-8:8)/2)),
            VC.resid = 1,
            nS = 2, nVc = 2, nSim = 10000, seed = 1, 
            neg.VC = TRUE, 
            row.MS = rowNames[c(8, 11)])

sim4Tags1$tempS = rbind(sim4Tags1$tempS, xx$tempS)
 
edfPlot(sim4Tags1)


sim4Tagstemp = sim4Tags1
 
sim4Tagstemp$tempS = sim4Tagstemp$tempS[,-5]

pdf("CRD232Tag4.pdf", width = 12)
edfPlot(sim4Tagstemp)
dev.off()

