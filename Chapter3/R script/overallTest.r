#13/08/2013 03:00:21 PM
library(optimTE)

#Even trt and compare between different techincal replicates
design.df = optCRD(nTrt = 4, bRep  = 2, tRep  = 2, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)

output = matrix(paste(design.df$Ani, design.df$Trt, sep = ""), ncol = 4, byrow = TRUE)
rownames(output) = 1:nrow(output)
output = rbind(c(114,115,116,117), output)

write.csv(file = "output.csv", output)

design.df$Tag = c(114,115,116,117)[design.df$Tag]
write.csv(file = "C:/Users/Kevin/Dropbox/OptimalDesigns/CRD/designCRD4224.csv", design.df,  row.names = FALSE)


design.df = optCRD(nTrt = 4, bRep  = 2, tRep  = 2, nPlot = 8, iter  = 1000)
design.summary.CRD(design.df, FALSE)

output = matrix(paste(design.df$Ani, design.df$Trt, sep = ""), ncol = 8, byrow = TRUE)
rownames(output) = 1:nrow(output)
output = rbind(c(113:119, 121), output)

write.csv(file = "output.csv", output)

design.df$Tag = c(113:119, 121)[design.df$Tag]
write.csv(file = "C:/Users/Kevin/Dropbox/OptimalDesigns/CRD/designCRD4228.csv", design.df,  row.names = FALSE)


design.df = optCRD(nTrt = 3, bRep  = 2, tRep  = 2, nPlot = 8, iter  = 1000)


##################################################################################################

design.df = optCRD(nTrt = 2, bRep  = 2, tRep  = 2, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)


design.df = optCRD(nTrt = 2, bRep  = 2, tRep  = 3, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 2, bRep  = 2, tRep  = 4, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)

#Odd trt and compare between different techincal replicates
design.df = optCRD(nTrt = 3, bRep  = 2, tRep  = 2, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 3, bRep  = 2, tRep  = 4, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)


#Even trt and compare between different techincal replicates
design.df = optCRD(nTrt = 4, bRep  = 2, tRep  = 2, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 4, bRep  = 2, tRep  = 3, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 4, bRep  = 2, tRep  = 4, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)


design.df = optCRD(nTrt = 5, bRep  = 4, tRep  = 2, nPlot = 8, iter  = 1000)
design.summary.CRD(design.df, FALSE)


design1.df = optCRD(nTrt = 5, bRep  = 8, tRep  = 2, nPlot = 8, iter  = 1000)
design.summary.CRD(design1.df, FALSE)

###########################################################################################

design.df = optCRD(nTrt = 8, bRep  = 6, tRep  = 2, nPlot = 96, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design1.df = optCRD(nTrt = 8, bRep  = 6, tRep  = 2, nPlot = 48, iter  = 1000)
design.summary.CRD(design1.df, FALSE)

design2.df = optCRD(nTrt = 8, bRep  = 6, tRep  = 2, nPlot = 24, iter  = 1000)
design.summary.CRD(design2.df, FALSE)

design3.df = optCRD(nTrt = 8, bRep  = 6, tRep  = 2, nPlot = 12, iter  = 1000)
design.summary.CRD(design3.df, FALSE)

design4.df = optCRD(nTrt = 8, bRep  = 6, tRep  = 2, nPlot = 6, iter  = 1000)
design.summary.CRD(design4.df, FALSE)



#Even trt and compare between different techincal replicates
design.df = optCRD(nTrt = 8, bRep  = 3, tRep  = 2, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 8, bRep  = 4, tRep  = 2, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 8, bRep  = 5, tRep  = 2, nPlot = 4, iter  = 10000)
design.summary.CRD(design.df, FALSE)


design.df = optCRD(nTrt = 8, bRep  = 10, tRep  = 2, nPlot = 4, iter  = 10000)
design.summary.CRD(design.df, FALSE)




#Even trt and compare between different techincal replicates
design.df = optCRD(nTrt = 4, bRep  = 3, tRep  = 2, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, latex = TRUE)


out <- matrix(paste0(design.df$Ani, design.df$Trt), ncol = 4, nrow = 6, byrow = TRUE)
xtable::xtable(out)

summaryAovTwoPhase(design.df, blk.str1 = "Ani", "Run", "Tag + Trt", decimal = TRUE)

summaryAovTwoPhase(design.df[-c(21:24),], blk.str1 = "Ani", "Run", "Tag + Trt", latex = TRUE)


summaryAovTwoPhase(design.df[-c(17:24),], blk.str1 = "Ani", "Run", "Tag + Trt", latex = TRUE)




design.df = optCRD(nTrt = 6, bRep  = 4, tRep  = 2, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df)

design.df = optCRD(nTrt = 6, bRep  = 5, tRep  = 2, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df)


design.df = optCRD(nTrt = 6, bRep  = 10, tRep  = 2, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df)




#Even trt and compare between different techincal replicates

design.df = optCRD(nTrt = 2, bRep  = 4, tRep  = 3, nPlot = 8, iter  = 1000)
design.summary.CRD(design.df, FALSE)



design.df = optCRD(nTrt = 2, bRep  = 4, tRep  = 4, nPlot = 8, iter  = 1000)

design.summary.CRD(design.df, FALSE)


#Odd trt and compare between different techincal replicates
design.df = optCRD(nTrt = 3, bRep  = 4, tRep  = 2, nPlot = 8, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 3, bRep  = 8, tRep  = 3, nPlot = 8, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 3, bRep  = 4, tRep  = 4, nPlot = 8, iter  = 1000)
design.summary.CRD(design.df, FALSE)

#Odd trt and compare between different techincal replicates
design.df = optCRD(nTrt = 6, bRep  = 2, tRep  = 4, nPlot = 8, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 6, bRep  = 2, tRep  = 4, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 6, bRep  = 2, tRep  = 3, nPlot = 4, iter  = 1000)
design.summary.CRD(design.df, FALSE)


design.df = optCRD(nTrt = 6, bRep  = 8, tRep  = 3, nPlot = 8, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 6, bRep  = 10, tRep  = 3, nPlot = 8, iter  = 1000)
design.summary.CRD(design.df, FALSE)

#Even trt and compare between different techincal replicates
design.df = optCRD(nTrt = 8, bRep  = 2, tRep  = 2, nPlot = 8, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 8, bRep  = 2, tRep  = 3, nPlot = 8, iter  = 1000)
design.summary.CRD(design.df, FALSE)

design.df = optCRD(nTrt = 8, bRep  = 2, tRep  = 4, nPlot = 8, iter  = 1000)
design.summary.CRD(design.df, FALSE)

###########################################################################################

optRBD.combine = function(nTrt, bRep, nCag, tRep, nPlot, iter, confoundCag = FALSE) {

    design.df = optCRD(nTrt = nTrt, bRep = bRep, tRep = tRep, nPlot = nPlot,
        iter = iter/10)
    n = nrow(design.df)
    nBlk = nlevels(design.df$Run)
    nPlot = nlevels(design.df$Tag)
    Run.mat = with(design.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
    Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
    Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))
    Pb = projMat(Run.mat)
    Pb1 = projMat(Tag.mat)

    (ave.eff = test.CRD(X.trt = Trt.mat, (diag(n) - Pb) %*% (diag(n) - Pb1),
        Rep = n/ncol(Trt.mat))$ave.eff)

    return(optRBD(nTrt = nTrt, bRep = bRep, nCag = nCag, tRep = tRep, nPlot = nPlot,
        iter = iter, upperValue = ave.eff, confoundCag = confoundCag))
} 
###########################################################################################

design.df = optRBD.combine(nTrt = 3, bRep = 4, nCag = 2, tRep = 2, nPlot = 4, iter = 1000)

design.summary.RBD(design.df)

design1.df = optRBD.combine(nTrt = 3, bRep = 4, nCag = 2, tRep = 2, nPlot = 4, iter = 1000, 
                            confoundCag = TRUE)

design.summary.RBD(design1.df)

design.df = optRBD(nTrt = 3, bRep = 4, nCag = 2, tRep = 2, nPlot = 4,
       iter = 100, confoundCag = FALSE)

design.summary.RBD(design.df)

design1.df = optRBD(nTrt = 3, bRep = 4, nCag = 2, tRep = 2, nPlot = 4,
                   iter = 100, confoundCag = TRUE)
design.summary.RBD(design1.df)



design.df = optRBD(nTrt = 4, bRep = 4, nCag = 4, tRep = 2, nPlot = 4,
                   iter = 100, confoundCag = FALSE)

design.summary.RBD(design.df, latex = TRUE)


design1.df = optRBD(nTrt = 4, bRep = 2, nCag = 2, tRep = 2, nPlot = 4,
                    iter = 100, confoundCag = TRUE)
design.summary.RBD(design1.df, latex = TRUE)


design.df = optRBD(nTrt = 4, bRep = 2, nCag = 2, tRep = 2, nPlot = 4,
                   iter = 100, confoundCag = FALSE)

design.summary.RBD(design.df)


design1.df = optRBD(nTrt = 4, bRep = 4, nCag = 4, tRep = 2, nPlot = 4,
                    iter = 100, confoundCag = TRUE)
design.summary.RBD(design1.df)

design.df = optRBD.combine(nTrt = 8, bRep = 6, nCag = 6, tRep = 2, nPlot = 4, iter = 1000)

design.summary.RBD(design.df)

output = matrix(paste(design.df$Cag, design.df$Ani, design.df$Trt, sep = ""), 
	ncol = 4, byrow = TRUE)
	
rownames(output) = 1:nrow(output)
output = rbind(c(114,115,116,117), output)

write.csv(file = "output.csv", output)

design.df$Tag = c(114,115,116,117)[design.df$Tag]
write.csv(file = "C:/Users/Kevin/Dropbox/OptimalDesigns/RCBD/designRCBD86624.csv", design.df,  row.names = FALSE)

#To compare
design1.df = optRBD.combine(nTrt = 8, bRep =6, nCag = 6, tRep = 2, nPlot = 4, iter = 3000, confoundCag = TRUE)

design.summary.RBD(design1.df)

output = matrix(paste(design1.df$Cag, design1.df$Ani, design1.df$Trt, sep = ""), 
	ncol = 4, byrow = TRUE)
	
rownames(output) = 1:nrow(output)
output = rbind(c(114,115,116,117), output)

write.csv(file = "output.csv", output)

design1.df$Tag = c(114,115,116,117)[design1.df$Tag]
write.csv(file = "C:/Users/Kevin/Dropbox/OptimalDesigns/RCBD/designRCBD86624EDF.csv", design1.df, row.names = FALSE)


#####################################################

design.df = optRBD.combine(nTrt = 8, bRep = 4, nCag = 4, tRep = 2, nPlot = 8, iter = 1000)
design.summary.RBD(design.df)

write.csv(file = "designRCBD84428.csv", design.df)

design1.df = optRBD.combine(nTrt = 8, bRep = 4, nCag = 4, tRep = 2, nPlot = 8, iter = 1000, 
                            confoundCag = TRUE)

write.csv(file = "designRCBD84428EDF.csv", design1.df)

design.summary.RBD(design1.df)


design.df = optRBD(nTrt = 8, bRep = 4, nCag = 4, tRep = 2, nPlot = 8, iter = 1000)
design.summary.RBD(design.df)



design.df = optRBD.combine(nTrt = 6, bRep = 2, nCag = 2, tRep = 2, nPlot = 4, iter = 1000)
design.summary.RBD(design.df)

design1.df = optRBD.combine(nTrt = 6, bRep = 2, nCag = 2, tRep = 2, nPlot = 4, 
                           iter = 1000, confoundCag = TRUE)
design.summary.RBD(design1.df)




design.df = optRBD.combine(nTrt = 6, bRep =8, nCag = 8, tRep = 2, nPlot = 4, iter = 1000)

design.summary.RBD(design.df)

output = matrix(paste(design.df$Cag, design.df$Ani, design.df$Trt, sep = ""), 
	ncol = 4, byrow = TRUE)
	
rownames(output) = 1:nrow(output)
output = rbind(c(114,115,116,117), output)

write.csv(file = "output.csv", output)

design.df$Tag = c(114,115,116,117)[design.df$Tag]
write.csv(file = "C:/Users/kcha193/Dropbox/OptimalDesigns/RCBD/designRCBD68824.csv", design.df,  row.names = FALSE)


design1.df = optRBD.combine(nTrt = 6, bRep = 8, nCag = 8, tRep = 2, nPlot = 4, iter = 1000, confoundCag = TRUE)

design.summary.RBD(design1.df)

output = matrix(paste(design1.df$Cag, design1.df$Ani, design1.df$Trt, sep = ""), 
	ncol = 4, byrow = TRUE)
	
rownames(output) = 1:nrow(output)
output = rbind(c(114,115,116,117), output)

write.csv(file = "output.csv", output)

design1.df$Tag = c(114,115,116,117)[design1.df$Tag]
write.csv(file = "C:/Users/kcha193/Dropbox/OptimalDesigns/RCBD/designRCBD68824EDF.csv", design1.df, row.names = FALSE)


####
#8-plex!!


design.df = optRBD(nTrt = 6, bRep = 6, nCag = 3, tRep = 2, nPlot = 8, iter = 1000)

design.summary.RBD(design.df, FALSE)


output = matrix(paste(design.df$Cag, design.df$Ani, design.df$Trt, sep = ""), ncol = 8, byrow = TRUE)
rownames(output) = 1:nrow(output)
output = rbind(c(113:119, 121), output)

write.csv(file = "output.csv", output)

design.df$Tag = c(113:119, 121)[design.df$Tag]
write.csv(file = "C:/Users/Kevin/Dropbox/OptimalDesigns/RCBD/designRCBD66328.csv", design.df,  row.names = FALSE)


#To compare
design1.df = optRBD.combine(nTrt = 6, bRep = 6, nCag = 6, tRep = 2, nPlot = 8, iter = 10000, confoundCag = TRUE)

design.summary.RBD(design1.df, FALSE)


output = matrix(paste(design1.df$Cag, design1.df$Ani, design1.df$Trt, sep = ""), ncol = 8, byrow = TRUE)
rownames(output) = 1:nrow(output)
output = rbind(c(113:119, 121), output)

write.csv(file = "output.csv", output)

design1.df$Tag = c(113:119, 121)[design1.df$Tag]
write.csv(file = "C:/Users/kcha193/Dropbox/OptimalDesigns/RCBD/designRCBD66628EDF.csv", design1.df,  row.names = FALSE)

######################################################################################
design1.df = optRBD.combine(nTrt = 6, bRep = 6, nCag = 2, tRep = 2, nPlot = 8, iter = 1000, confoundCag = TRUE)

design.summary.RBD(design1.df, FALSE)


output = matrix(paste(design1.df$Cag, design1.df$Ani, design1.df$Trt, sep = ""), ncol = 8, byrow = TRUE)
rownames(output) = 1:nrow(output)
output = rbind(c(113:119, 121), output)

write.csv(file = "output.csv", output)

design1.df$Tag = c(113:119, 121)[design1.df$Tag]
write.csv(file = "C:/Users/Kevin/Dropbox/OptimalDesigns/RCBD/designRCBD66228EDF.csv", design1.df,  row.names = FALSE)



###########################################################################################
#Difficult ones !!!


design.df = optRBD.combine(nTrt = 2, bRep = 6, nCag = 6, tRep = 2, nPlot = 4, iter = 1000)
design.summary.RBD(design.df)


output = matrix(paste(design.df$Cag, design.df$Ani, design.df$Trt, sep = ""), 
	ncol = 4, byrow = TRUE)
	
rownames(output) = 1:nrow(output)
output = rbind(c(114,115,116,117), output)

write.csv(file = "output.csv", output)

design.df$Tag = c(114,115,116,117)[design.df$Tag]
write.csv(file = "C:/Users/Kevin/Dropbox/OptimalDesigns/RCBD/designRCBD26624.csv", design.df,  row.names = FALSE)




design.df = optRBD.combine(nTrt = 2, bRep = 4, nCag = 2, tRep = 2, nPlot = 8, iter = 1000)
design.summary.RBD(design.df)


design.summary.RBD(design.df, FALSE)


output = matrix(paste(design.df$Cag, design.df$Ani, design.df$Trt, sep = ""), ncol = 8, byrow = TRUE)
rownames(output) = 1:nrow(output)
output = rbind(c(113:119, 121), output)

write.csv(file = "output.csv", output)

design.df$Tag = c(113:119, 121)[design.df$Tag]
write.csv(file = "C:/Users/Kevin/Dropbox/OptimalDesigns/BIBD/designRCBD24228.csv", design.df,  row.names = FALSE)

######################################################################################################





design1.df = optRBD.combine(nTrt = 2, bRep = 4, nCag = 4, tRep = 2, nPlot = 4, iter = 1000)
design.summary.RBD(design1.df)



design1.df = optRBD.combine(nTrt = 2, bRep = 4, nCag = 4, tRep = 3, nPlot = 4, iter = 1000)
design.summary.RBD(design1.df)



design1.df = optRBD.combine(nTrt = 3, bRep = 2, nCag = 2, tRep = 2, nPlot = 4, iter = 1000)
design.summary.RBD(design1.df)



design1.df = optRBD.combine(nTrt = 6, bRep = 4, nCag = 4, tRep = 2, nPlot = 4, iter = 1000)
design.summary.RBD(design1.df)

design2.df = optRBD.combine(nTrt = 6, bRep = 4, nCag = 4, tRep = 2, nPlot = 4, iter = 1000, confoundCag = TRUE)
design.summary.RBD(design2.df)


design1.df = optRBD.combine(nTrt = 6, bRep = 4, nCag = 2, tRep = 2, nPlot = 4, iter = 1000)
design.summary.RBD(design1.df)

design2.df = optRBD.combine(nTrt = 6, bRep = 4, nCag = 2, tRep = 2, nPlot = 4, iter = 1000, confoundCag = TRUE)
design.summary.RBD(design2.df)


design1.df = optRBD.combine(nTrt = 6, bRep = 8, nCag = 8, tRep = 2, nPlot = 4, iter = 10000)

design2.df = optRBD.combine(nTrt = 6, bRep = 8, nCag = 8, tRep = 2, nPlot = 4, iter = 10000, confoundCag = TRUE)

design.summary.RBD(design1.df, FALSE)
design.summary.RBD(design2.df, FALSE)


design1.df = optRBD.combine(nTrt = 6, bRep = 6, nCag = 6, tRep = 2, nPlot = 8, iter = 1000)

design2.df = optRBD.combine(nTrt = 6, bRep = 6, nCag = 6, tRep = 2, nPlot = 8, iter = 1000, confoundCag = TRUE)

design.summary.RBD(design1.df, FALSE)
design.summary.RBD(design2.df, FALSE)

design1.df = optRBD.combine(nTrt = 6, bRep = 6, nCag = 3, tRep = 2, nPlot = 8, iter = 1000)

design.summary.RBD(design1.df, FALSE)

design1.df = optRBD.combine(nTrt = 6, bRep = 10, nCag = 10, tRep = 2, nPlot = 8, iter = 10000)

design2.df = optRBD.combine(nTrt = 6, bRep = 10, nCag = 10, tRep = 2, nPlot = 8, iter = 10000, confoundCag = TRUE)

design.summary.RBD(design1.df, FALSE)
design.summary.RBD(design2.df, FALSE)



design1.df = optRBD.combine(nTrt = 8, bRep = 10, nCag = 10, tRep = 2, nPlot = 8, iter = 1000)

design2.df = optRBD.combine(nTrt = 8, bRep = 10, nCag = 10, tRep = 2, nPlot = 8, iter = 1000, confoundCag = TRUE)

design.summary.RBD(design1.df, FALSE)
design.summary.RBD(design2.df, FALSE)


design1.df = optRBD.combine(nTrt = 8, bRep = 8, nCag = 8, tRep = 2, nPlot = 4, iter = 1000)

design2.df = optRBD.combine(nTrt = 8, bRep = 8, nCag = 8, tRep = 2, nPlot = 4, iter = 1000, confoundCag = TRUE)

design.summary.RBD(design1.df, FALSE)
design.summary.RBD(design2.df, FALSE)


#################################################################################################

design.df = optBIBD(nTrt = 8, bRep = 7, nCag = 8, tRep = 2, nPlot = 4, resDF = 29, iter = 3000, confoundCag = TRUE)

design.summary.RBD(design.df, FALSE)


output = matrix(paste(design.df$Cag, design.df$Ani, design.df$Trt, sep = ""), 
	ncol = 4, byrow = TRUE)
	
rownames(output) = 1:nrow(output)
output = rbind(c(114,115,116,117), output)

write.csv(file = "output.csv", output)

design.df$Tag = c(114,115,116,117)[design.df$Tag]
write.csv(file = "C:/Users/Kevin/Dropbox/OptimalDesigns/BIBD/designBIBD87824.csv", design.df,  row.names = FALSE)




design.df = optBIBD(nTrt = 7, bRep = 6, nCag = 7, tRep = 2, nPlot = 4, iter = 10000, confoundCag = FALSE)

design.summary.RBD(design.df, FALSE)


output = matrix(paste(design.df$Cag, design.df$Ani, design.df$Trt, sep = ""), 
	ncol = 4, byrow = TRUE)
	
rownames(output) = 1:nrow(output)
output = rbind(c(114,115,116,117), output)

write.csv(file = "output.csv", output)

design.df$Tag = c(114,115,116,117)[design.df$Tag]
write.csv(file = "C:/Users/Kevin/Dropbox/OptimalDesigns/BIBD/designBIBD76724.csv", design.df,  row.names = FALSE)

write.csv(file = "C:/Users/kcha193/Dropbox/OptimalDesigns/BIBD/output.csv", output)

write.csv(file = "C:/Users/kcha193/Dropbox/OptimalDesigns/BIBD/designBIBD76724.csv", design.df,  row.names = FALSE)

$can.eff
[1] 0.9209592 0.8736100 0.8628041 0.8398796 0.8138405 0.7802297

0.921 0.874 0.863 0.840 0.814 0.780

$ave.eff
[1] 0.8461938


	
	
design.df = optBIBD(nTrt = 8, bRep = 7, nCag = 8, tRep = 2, nPlot = 8, iter = 3000, confoundCag = TRUE)

design.summary.RBD(design.df, FALSE)

output = matrix(paste(design.df$Cag, design.df$Ani, design.df$Trt, sep = ""), ncol = 8, byrow = TRUE)
rownames(output) = 1:nrow(output)
output = rbind(c(113:119, 121), output)

write.csv(file = "output.csv", output)

design.df$Tag = c(113:119, 121)[design.df$Tag]
write.csv(file = "C:/Users/Kevin/Dropbox/OptimalDesigns/BIBD/designBIBD87828.csv", design.df,  row.names = FALSE)


#################################################################################################
design.df = optBIBD(nTrt = 4, bRep = 3, nCag = 4, tRep = 2, nPlot = 4,
iter = 1000, confoundCag = TRUE)

design.summary.RBD(design.df, FALSE)
                 
design.df = optBIBD(nTrt = 4, bRep = 3, nCag = 4, tRep = 2, nPlot = 8,
iter = 1000, confoundCag = TRUE)

design.summary.RBD(design.df, FALSE)

design.df = optBIBD(nTrt = 5, bRep = 4, nCag = 5, tRep = 2, nPlot = 4,
iter = 1000)

design.summary.RBD(design.df, FALSE)

design1.df = optRBD.combine(nTrt = 5, bRep = 4, nCag = 4, tRep = 2, nPlot = 4,
iter = 1000, confoundCag = TRUE)

design.summary.RBD(design1.df, FALSE)

                 
design.df = optBIBD(nTrt = 5, bRep = 4, nCag = 5, tRep = 2, nPlot = 8,
iter = 1000)

design.summary.RBD(design.df, FALSE)



design.df = optBIBD(nTrt = 6, bRep = 5, nCag = 6, tRep = 2, nPlot = 4,
iter = 3000, confoundCag = TRUE)

design.summary.RBD(design.df, FALSE)
            

speDes = matrix(c("a", "b", "c",
                  "a", "d", "e",
                  "a", "f", "g",
                  "b", "d", "f",
                  "b", "e", "g",
                  "c", "d", "g",
                  "c", "e", "f"), nrow = 7, ncol = 3, byrow = TRUE)

				  
design.df = optBIBD(nTrt = 7, bRep = 4, nCag = 7, tRep = 2, nPlot = 8, iter = 1000,
    confoundCag = FALSE, speDes = speDes)
			

design.df = optBIBD(nTrt = 7, bRep = 6, nCag = 7, tRep = 2, nPlot = 4,
iter = 100000)

design.summary.RBD(design.df, FALSE)


design.df = optBIBD(nTrt = 8, bRep = 7, nCag = 8, tRep = 2, nPlot = 4,
iter = 3000, confoundCag = TRUE)

design.summary.RBD(design.df, FALSE)


design.df = optBIBD(nTrt = 8, bRep = 7, nCag = 8, tRep = 2, nPlot = 8,
iter = 1000, confoundCag = TRUE)

design.summary.RBD(design.df, FALSE)


