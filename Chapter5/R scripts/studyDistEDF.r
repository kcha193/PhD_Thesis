


summaryEDF = cbind(apply(REAL.EDF[-1,], 2, function(x) median(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) mean(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) mean(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(REML.EDF[-1,], 2, function(x) median(x[which(!is.infinite(x))], na.rm = TRUE)),
apply(LC.EDF[-1,], 2, function(x) median(x[which(!is.infinite(x))], na.rm = TRUE)))

colnames(summaryEDF) = c("True.EDF", "REML.EDF.mean", "LC.EDF.mean", "REML.EDF.median", "LC.EDF.median")

summaryEDF

summaryVC = cbind(apply(REAL.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) mean(x, na.rm = TRUE)),
apply(REML.VC[-1,], 2, function(x) median(x, na.rm = TRUE)),
apply(LC.VC[-1,], 2, function(x) median(x, na.rm = TRUE)))

colnames(summaryVC) = c("True.EDF", "REML.VC.mean", "LC.VC.mean",
                                     "REML.VC.median", "LC.VC.median")
summaryVC


histEDF =
function(ind.REML, ind.LC, bin = 0.1, main) {
  data.EDF = data.frame(Method = rep(c("REML", "LC"), c(length(REML.EDF[ind.REML,3]), length(LC.EDF[ind.LC,3]))),
                        EDF = c(REML.EDF[ind.REML,3], LC.EDF[ind.LC,3])  )

  vline.mean <- data.frame(Method = c("REML", "LC"), Mean = c(mean(REML.EDF[ind.REML,3], na.rm = TRUE), mean(LC.EDF[ind.LC,3], na.rm = TRUE)))
  vline.median <- data.frame(Method = c("REML", "LC"), Median = c(median(REML.EDF[ind.REML,3], na.rm = TRUE), median(LC.EDF[ind.LC,3], na.rm = TRUE)))

  print(cbind(vline.mean, Median = vline.median[,-1]))

  pg = ggplot(data.EDF, aes(x =EDF)) +
      ggtitle(main) +
      geom_histogram( binwidth = bin)+
      facet_wrap(~ Method) +
      geom_vline(aes(xintercept = mean(REAL.EDF[-1,3])), col = "Red") +
      geom_vline(aes(xintercept = Mean), col = "Darkgreen", vline.mean)+
      geom_vline(aes(xintercept = Median), col = "Blue", vline.median)

  pg
}

bin = VC.resid

grid.arrange( histEDF (-1, -1, bin = bin, main = "Original"),
              histEDF(-c(1,which(REML.VC[,2]<0)), -c(1,which(LC.VC[,2]<0)), bin = bin, main = "Negative Ani VC Removed"),
              histEDF(-c(1,which(REML.VC[,3]<0)), -c(1,which(LC.VC[,3]<0)), bin = bin, main = "Negative Run VC Removed"),
              histEDF( -c(1,union(which(REML.VC[,2]<0), which(REML.VC[,3]<0))),
              -c(1,union(which(LC.VC[,2]<0),which(LC.VC[,3]<0))), bin = bin, main = "Negative Ani and Run VC Removed"), ncol=1)


bin = 0.1


grid.arrange( histEDF (-1, -1, bin = bin, main = "Original"),
              histEDF(-c(1,which(REML.VC[,2]==0)), -c(1,which(LC.VC[,2]==0)), bin = bin, main = "Zero Ani VC Removed"),
              histEDF(-c(1,which(REML.VC[,3]==0)), -c(1,which(LC.VC[,3]==0)), bin = bin, main = "Zero Run VC Removed"),
              histEDF( -c(1,union(which(REML.VC[,2]==0), which(REML.VC[,3]==0))),
              -c(1,union(which(LC.VC[,2]==0),which(LC.VC[,3]==0))), bin = bin, main = "Zero Ani and Run VC Removed"), ncol=1)

