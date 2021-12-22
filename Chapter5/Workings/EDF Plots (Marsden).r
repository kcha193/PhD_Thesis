## EDF profiles stored in 3 different dataframes.
## The structure of these is:

        |     sa2/s2
sp2/s2  | <--- (-8:8)/2 --->
--------+--------------------
  0     |
  0.01  |
  0.25  |
  1     |
  4     |
100     |
Fixed   |

## Note that there is no Fixed row for ALD using J&W method

## Want to put these into a single data frame
## Need 3 columns: EDF, sp2/s2 value (factor), sa2/2 value (factor)


## Put these into one dataframe.

edf.df <- local(
{
   edf <- rbind(stack(data.frame(t(outpmds[-c(2:4),]))), # Eliminate case sa2/s2=0.01
                stack(data.frame(t(outp[-c(2:4),]))),
                stack(data.frame(t(outplc[-c(2:4),]))))
   edf.mat <- as.matrix(edf)
   edf$spR <- rep(c(10^((-8:8)/2)), times=11)
   edf$saR <- gsub("X","",edf.mat[,"ind"])
#   edf$saR <- factor(gsub("X",expression(paste(sigma^2[alpha]," / ", sigma^2, " = "), edf.mat[,"ind"]),
   edf$saR <- factor(edf$saR, levels=unique(edf$saR)[c(7,1:6)],ordered=T)
   edf <- edf[,-2]
   dimnames(edf)[[2]][1] <- "EDF"
   edf$Method <- factor(rep(1:3,times=c(68,68,51)),
                       labels=c("MDS/REML","ALD/REML", "ALD/LC"), ordered=T)
   edf
})

## Code for Marsden proposal

require(lattice)
trellis.device(theme=col.whitebg())

trellis.par.set(layout.heights = list(strip = 1.25))
xyplot(EDF~spR|saR, group=Method, data=edf.df, lty=1:3, type="l", #as.table=T,
   col=c("black", "red", "dark blue"), lwd=2,
   scales=list(x=list(at=10^(-3:3), labels=c(0.001,.01,.1,1,10,100,1000), log=TRUE)),
   xlab=expression(paste(sigma[pi]^2," / ",sigma^2)),
   key =
      list(text = list(paste(c("Design/Method:  ","", ""), levels(edf.df$Method))),
           space = "top", 
           lines = list(lty=1:3, col=c("black", "red", "dark blue"), lwd=2),
           columns=3, rep=FALSE
           ),
   index.cond=list(c(3:4,1:2)),
   strip = function(factor.levels, par.strip.text, ...)
             {strip.default(factor.levels = expression(paste(sigma[alpha]^2, " / ", sigma^2, " = 0"),
                                                    paste(sigma[alpha]^2, " / ", sigma^2, " = 4"),
                                                    paste(sigma[alpha]^2, " / ", sigma^2, " = 100"),
                                                    "Arrays Fixed"),
     par.strip.text = trellis.par.get("layout.heights"), ...)}
)
dev.print(device=postscript,width = 10 * 140, height = 6.5 * 140)


## OLD code for Revision 1

trellis.par.set(layout.heights = list(strip = 1.25))
xyplot(EDF~spR|saR, group=Method, data=edf.df, lty=1:3, type="l", as.table=T,
   col=c("black", "red", "dark blue"), lwd=2,
   scales=list(x=list(at=10^(-3:3), labels=c(0.001,.01,.1,1,10,100,1000), log=TRUE)),
   xlab=expression(paste(sigma[pi]^2," / ",sigma^2)), 
   key =
      list(text = list(paste(c("Design/Method:  ","", ""), levels(edf.df$Method))),
           space = "top", 
           lines = list(lty=1:3, col=c("black", "red", "dark blue"), lwd=2),
           columns=3, rep=FALSE
           ),
   index.cond=list(c(2:7,1)),
   strip = strip.custom(strip.names=TRUE,var.name=expression(paste(sigma[alpha]^2, " / ", sigma^2)), sep=" = ",
     par.strip.text = trellis.par.get("layout.heights"))
)
trellis.par.set(layout.heights = list(strip = 1))


xyplot(EDF~spR|saR, group=Method, data=edf.df, lty=1:3, type="l")



trellis.par.set(layout.heights = list(strip = 1.25))
xyplot(EDF~spR|saR, group=Method, data=edf.df, lty=1:3, type="l", #as.table=T,
   col=c("black", "red", "dark blue"), lwd=2,
   scales=list(x=list(at=10^(-3:3), labels=c(0.001,.01,.1,1,10,100,1000), log=TRUE)),
   xlab=expression(paste(sigma[pi]^2," / ",sigma^2)), 
   key =
      list(text = list(paste(c("Design/Method:  ","", ""), levels(edf.df$Method))),
           space = "top", 
           lines = list(lty=1:3, col=c("black", "red", "dark blue"), lwd=2),
           columns=3, rep=FALSE,
           ),
   index.cond=list(c(5:6,1,2:4)),
   strip = function(factor.levels, par.strip.text, ...)
             {strip.default(factor.levels = expression("Arrays Fixed",
                                                    paste(sigma[alpha]^2, " / ", sigma^2, " = 0"),
                                                    paste(sigma[alpha]^2, " / ", sigma^2, " = 0.25"),
                                                    paste(sigma[alpha]^2, " / ", sigma^2, " = 1"),
                                                    paste(sigma[alpha]^2, " / ", sigma^2, " = 4"),
                                                    paste(sigma[alpha]^2, " / ", sigma^2, " = 100")),
     par.strip.text = trellis.par.get("layout.heights"), ...)}
)
dev.print(device=postscript,width = 10 * 140, height = 6.5 * 140)
