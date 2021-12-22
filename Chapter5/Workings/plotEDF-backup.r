###########################################################################################
#Plot

  tempS.bwAni = sim$tempS[grep("WithinBetweenAniResidual", rownames(sim$tempS)),]

  mm = melt(tempS.bwAni)

  new.temp = mm[which(!is.na(match(mm$variable, c("REML.EDF", "LC.EDF", "REAL.EDF")))),]

  new.temp$value = round(new.temp$value, digits = 4)
  colnames(new.temp)[4] = "EDF"

   trellis.par.set(layout.heights = list(strip = 1.25),
          superpose.line = list(col=c("black", "red", "blue"), lwd=2))
   xyplot(EDF ~ gamma.ani|gamma.run, group = variable, type="l",
   col=c("black", "red", "blue"), lwd=2,
    xlab=expression(paste(sigma[A]^2," / ",sigma^2)),
      key  =
      simpleKey(text = paste(c("Method:  ","", ""), levels(new.temp$variable)),
          points = FALSE, lines = TRUE, cex = 1, columns=3),
    strip = function(factor.levels, par.strip.text, ...)
             {strip.default(factor.levels = expression(
                             paste(sigma[R]^2, " / ", sigma^2, " = 1e-4"),
                             paste(sigma[R]^2, " / ", sigma^2, " = 1e-3"),
                             paste(sigma[R]^2, " / ", sigma^2, " = 1e-2"),
                             paste(sigma[R]^2, " / ", sigma^2, " = 0.01"),
                             paste(sigma[R]^2, " / ", sigma^2, " = 1"),
                             paste(sigma[R]^2, " / ", sigma^2, " = 10"),
                             paste(sigma[R]^2, " / ", sigma^2, " = 100"),
                             paste(sigma[R]^2, " / ", sigma^2, " = 1e3"),
                             paste(sigma[R]^2, " / ", sigma^2, " = 1e4")),
     par.strip.text = trellis.par.get("layout.heights"), ...)},
   data = new.temp )

###########################################################################################

    tempS.bwRun = sim$tempS[grep("Between Run", rownames(sim$tempS)),]


      mm = melt(tempS.bwRun)

  new.temp = mm[which(!is.na(match(mm$variable, c("REML.EDF", "LC.EDF", "REAL.EDF")))),]

  new.temp$value = round(new.temp$value, digits = 4)

  colnames(new.temp)[4] = "EDF"

   trellis.par.set(layout.heights = list(strip = 1.25),
          superpose.line = list(col=c("black", "red", "blue"), lwd=2))
   xyplot(EDF ~ gamma.run|gamma.ani, group = variable, type="l",
   col=c("black", "red", "blue"), lwd=2,
    xlab=expression(paste(sigma[R]^2," / ",sigma^2)),
      key  =
      simpleKey(text = paste(c("Method:  ","", ""), levels(new.temp$variable)),
          points = FALSE,
           lines = TRUE,
           cex = 1,
           columns=3, #rep=FALSE,
           ),
    strip = function(factor.levels, par.strip.text, ...)
             {strip.default(factor.levels = expression(
                            paste(sigma[A]^2, " / ", sigma^2, " = 1e-4"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 1e-3"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 1e-2"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 0.01"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 1"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 10"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 100"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 1e3"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 1e4")),
     par.strip.text = trellis.par.get("layout.heights"), ...)},
   data = new.temp )

   ###########################################################################################

    tempS.resid = sim$tempS[grep("WithinResidualResidual", rownames(sim$tempS)),]


      mm = melt(tempS.resid)

  new.temp = mm[which(!is.na(match(mm$variable, c("REML.EDF", "LC.EDF", "REAL.EDF")))),]

  new.temp$value = round(new.temp$value, digits = 4)

  colnames(new.temp)[4] = "EDF"

   trellis.par.set(layout.heights = list(strip = 1.25),
          superpose.line = list(col=c("black", "red", "blue"), lwd=2))
   xyplot(EDF ~ gamma.run|gamma.ani, group = variable, type="l",
   col=c("black", "red", "blue"), lwd=2,
    xlab=expression(paste(sigma[R]^2," / ",sigma^2)),
      key  =
      simpleKey(text = paste(c("Method:  ","", ""), levels(new.temp$variable)),
          points = FALSE,
           lines = TRUE,
           cex = 1,
           columns=3, #rep=FALSE,
           ),
    strip = function(factor.levels, par.strip.text, ...)
             {strip.default(factor.levels = expression(
                            paste(sigma[A]^2, " / ", sigma^2, " = 1e-4"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 1e-3"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 1e-2"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 0.01"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 1"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 10"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 100"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 1e3"),
                             paste(sigma[A]^2, " / ", sigma^2, " = 1e4")),
     par.strip.text = trellis.par.get("layout.heights"), ...)},
   data = new.temp )




    tempV.ani = sim$tempV[grep("Ani", rownames(sim$tempV)),]

    tempV.run = tempV[grep("Run", rownames(tempV)),]

    mm = melt(tempV.ani)

    new.temp = mm[which(!is.na(match(mm$variable, c("REML.VC", "LC.VC", "REAL.VC")))),]

    new.temp$value = round(new.temp$value, digits = 4)


    xyplot(value ~ gamma.ani|gamma.run, group = variable, type="l",
          col=c("black", "red", "dark blue"), lwd=2,
          xlab=expression(paste(sigma[A]^2," / ",sigma^2)),

          strip = function(factor.levels, par.strip.text, ...)
             {strip.default(factor.levels = expression("Arrays Fixed",
                                  paste(sigma[R]^2, " / ", sigma^2, " = 0"),
                                  paste(sigma[R]^2, " / ", sigma^2, " = 0.25"),
                                  paste(sigma[R]^2, " / ", sigma^2, " = 1"),
                                  paste(sigma[R]^2, " / ", sigma^2, " = 4"),
                                  paste(sigma[R]^2, " / ", sigma^2, " = 100"),
                                  paste(sigma[R]^2, " / ", sigma^2, " = 1000")),
     par.strip.text = trellis.par.get("layout.heights"), ...)},
   data = new.temp )


}



s = 1
sp = c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 1, 10, 100, 1000, 10000, 100000)
sa = 0
r = 4

a = r*(s + 2*sp)^2
b = r * s^2 + (r - 1) * (s + 2 * sa)^2 + r*(s + 2*sa + 2*sp)^2
EDF = (r-1) * (1+ a/b)
EDF

plot(1:length(sp), EDF, type = "l", xaxt = "n")

lines(1:length(sp), EDF, col = "red")
lines(1:length(sp), EDF, col = "orange")
lines(1:length(sp), EDF, col = "yellow")
lines(1:length(sp), EDF, col = "green")
lines(1:length(sp), EDF, col = "blue")



axis(1, 1:length(sp), as.character(sp))





