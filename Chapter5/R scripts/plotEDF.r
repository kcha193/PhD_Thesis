###########################################################################################
#Plot

plotEDF =
function(sim, rowNames = "Within RunBetweenAniResidual", ...){
  require(reshape)
  require(lattice)
  
  tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)),]

  mm = melt(tempS.bwAni)

  new.temp = mm[which(!is.na(match(mm$variable, c("REML.EDF", "LC.EDF", "REAL.EDF")))),]

  new.temp$value = round(new.temp$value, digits = 4)
  colnames(new.temp)[4] = "EDF"


   suppressWarnings(trellis.par.set(layout.heights = list(strip = 1.25),
          superpose.line = list(col=c("black", "red", "blue"), lwd=2,lty=1:3)))     
                                                                                      
   interval = as.numeric(levels(new.temp$gamma.ani)[seq(1,nlevels(new.temp$gamma.ani),2)])
   
   new.temp$gamma.ani = as.numeric(as.character(new.temp$gamma.ani))
   
   xyplot(EDF ~ gamma.ani|gamma.run, group = variable, data = new.temp,lty=1:3, type="l",
   col=c("black", "red", "blue"), lwd=2,
   scales=list(x=list(at=interval, labels=interval, log=TRUE)),
    xlab=expression(paste(sigma[A]^2," / ",sigma^2)),
      key  =
      simpleKey(text = paste(c("Method:  ","", ""), levels(new.temp$variable)),
          points = FALSE, lines = TRUE, cex = 1, columns=3),
     strip = function(factor.levels, par.strip.text, ...)
    {strip.default(factor.levels =  
    sapply(levels(new.temp$gamma.r), 
    function(x) as.expression(substitute(list(sigma[R]^2/sigma^2==x), list(x=x) ))),
  par.strip.text = trellis.par.get("layout.heights"), ...)},   ... )
   

}

edfPlot =
function(sim, rowNames = "Within RunBetweenAniResidual", pTitle = "EDF plots"){
  
  
  tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)),]

  mm = melt(tempS.bwAni)

  names(mm)[3] = "Method"
  
  mm$gamma.ani =   as.numeric(as.character(mm$gamma.ani))
  #mm$gamma.run =  as.factor( paste("gamma.run =", mm$gamma.run))
  #levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))
  
  
   pg <- ggplot(mm, aes(x = log(gamma.ani), y=value), group = gamma.run)+ ylab("EDF") +
  xlab(expression(paste(sigma[A]^2," / ",sigma^2))) +
  ggtitle(pTitle) +
  scale_x_continuous( breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)],
                      labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) +
  geom_line(aes(colour = Method), size = 1, alpha=I(.5)) +
  facet_wrap(~  gamma.run ) +
  scale_colour_discrete(name = "Method")
     
 pg
}

edfPlotCompare =
function(sim, sim1, compare, method = "REML", rowNames = "Within RunBetweenAniResidual", pTitle = "EDF plots"){
  
  
  tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)), c(1,2, grep(method, colnames(sim$tempS)))]
  tempS.bwAni = cbind(tempS.bwAni, sim1$tempS[grep(rowNames, rownames(sim1$tempS)),grep(method, colnames(sim1$tempS))])
  
  colnames(tempS.bwAni) [c(3,4)]=  compare
  mm = melt(tempS.bwAni)  

  colnames(mm) [3] = "Design"
  
  mm$gamma.ani =   as.numeric(as.character(mm$gamma.ani))
  #mm$gamma.run =  as.factor( paste("gamma.run =", mm$gamma.run))
  #levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))

   pg <- ggplot(mm, aes(x = log(gamma.ani), y=value))+ ylab("EDF") +
  xlab(expression(paste(sigma[A]^2," / ",sigma^2))) +
  ggtitle(pTitle) +
  scale_x_continuous( breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)],
                      labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) +
  geom_line(aes(colour = Design), size = 1, alpha=I(.5)) +
  facet_wrap(~ gamma.run) +
  scale_colour_discrete(name = "Design")
     
 pg
}


VCPlot =
function(sim, rowNames = "e", pTitle = "VC plots"){
  
  
  tempS.bwAni = sim$tempV[grep(rowNames, rownames(sim$tempV)),]

  mm = melt(tempS.bwAni)

  names(mm)[3] = "Method"
  
  mm$gamma.ani =   as.numeric(as.character(mm$gamma.ani))
  #mm$gamma.run =  as.factor( paste("gamma.run =", mm$gamma.run))
  #levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))
  
  
   pg <- ggplot(mm, aes(x = log(gamma.ani), y=value), group = gamma.run)+ ylab("VC") +
  xlab(expression(paste(sigma[A]^2," / ",sigma^2))) +
  ggtitle(pTitle) +
  scale_x_continuous( breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)],
                      labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) +
  geom_line(aes(colour = Method), size = 1, alpha=I(.5)) +
  facet_wrap(~  gamma.run ) +
  scale_colour_discrete(name = "Method")
     
 pg
}

###########################################################################################


edfPlot1 = function(sim, rowNames = "Within RunBetweenCag:AniResidual", pTitle = "EDF plots") {
    
    tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)), ]
    
    mm = melt(tempS.bwAni)
    
    names(mm)[4] = "Method"
    
    mm$gamma.ani = as.numeric(as.character(mm$gamma.ani))
    # mm$gamma.run = as.factor( paste('gamma.run =', mm$gamma.run)) levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))
    
    level.run = levels(mm$gamma.run)
    
    pg = list()
    for (i in 1:length(level.run)) {
        temp = mm[which(mm$gamma.run == level.run[i]), ]
        
        pg[[i]] <- ggplot(temp, aes(x = log(gamma.ani), y = value)) + 
        ylab("EDF") + xlab(expression(paste(sigma[A]^2, " / ", sigma^2))) + 
        ggtitle(as.expression(substitute(list(sigma[B]^2/sigma^2 == x), list(x = level.run[i])))) + 
        scale_x_continuous(breaks = log(unique(temp$gamma.ani))[seq(1, 17, 2)], 
                labels = as.character(unique(round(temp$gamma.ani, 4)))[seq(1, 17, 2)]) + 
        geom_line(aes(colour = Method), size = 1, alpha = I(0.5)) + 
        facet_wrap(~gamma.cag) + 
        scale_colour_discrete(name = "Method")
    }
    
    pg
}


edfPlot1.comCag = function(sim, rowNames = "Within RunBetweenCag:AniResidual", method = "REML", pTitle = "EDF plots") {
    
    tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)), 
                    c("gamma.run", "gamma.cag", "gamma.ani", paste(method, ".EDF", sep = ""))]
    
    mm = melt(tempS.bwAni)
    
    #names(mm)[4] = "Method"
    
    mm$gamma.ani = as.numeric(as.character(mm$gamma.ani))
    # mm$gamma.run = as.factor( paste('gamma.run =', mm$gamma.run)) levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))
    
    if(any(is.na(mm$gamma.run))){
    level.temp =  levels(mm$gamma.run)
    temp =  as.character(mm$gamma.run)
    temp[which(is.na(temp))] = "Fixed"
    mm$gamma.run = factor(temp,  levels =  c(level.temp, "Fixed"))
   }
   
   if(any(is.na(mm$gamma.cag))){
    level.temp =  levels(mm$gamma.cag)
    temp =  as.character(mm$gamma.cag)
    temp[which(is.na(temp))] = "Fixed"
    mm$gamma.cag = factor(temp,  levels =  c(level.temp, "Fixed"))
   }
 
    pg <- ggplot(mm, aes(x = log(gamma.ani), y = value)) + 
    ylab("EDF") + xlab(expression(paste(sigma[A]^2, " / ", sigma^2))) + 
    ggtitle(pTitle) + 
    scale_x_continuous(breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)], 
            labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) + 
    geom_line(aes(colour = gamma.cag), size = 1, alpha = I(0.5)) + 
    facet_wrap(~gamma.run) + 
    scale_colour_discrete(name = "gamma.cag")
  
    
    pg
}
 
    
    
edfPlot1.comRun = function(sim, rowNames = "Within RunBetweenCag:AniResidual", method = "REML", pTitle = "EDF plots") {
    
    tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)), 
                    c("gamma.run", "gamma.cag", "gamma.ani", paste(method, ".EDF", sep = ""))]
    
    mm = melt(tempS.bwAni)
    
    #names(mm)[4] = "Method"
    
    mm$gamma.ani = as.numeric(as.character(mm$gamma.ani))
    # mm$gamma.run = as.factor( paste('gamma.run =', mm$gamma.run)) levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))
    
    if(any(is.na(mm$gamma.run))){
    level.temp =  levels(mm$gamma.run)
    temp =  as.character(mm$gamma.run)
    temp[which(is.na(temp))] = "Fixed"
    mm$gamma.run = factor(temp,  levels =  c(level.temp, "Fixed"))
   }
   
   if(any(is.na(mm$gamma.cag))){
    level.temp =  levels(mm$gamma.cag)
    temp =  as.character(mm$gamma.cag)
    temp[which(is.na(temp))] = "Fixed"
    mm$gamma.cag = factor(temp,  levels =  c(level.temp, "Fixed"))
   }
 
    pg <- ggplot(mm, aes(x = log(gamma.ani), y = value)) + 
    ylab("EDF") + xlab(expression(paste(sigma[A]^2, " / ", sigma^2))) + 
    ggtitle(pTitle) + 
    scale_x_continuous(breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)], 
            labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) + 
    geom_line(aes(colour = gamma.run), size = 1, alpha = I(0.5)) + 
    facet_wrap(~gamma.cag) + 
    scale_colour_discrete(name = "gamma.run")
  
    
    pg
}

plotEDF1 =
function(sim, rowNames = "Within RunBetweenCag:AniResidual", ...){
  require(reshape)
  require(lattice)
  
  tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)),]

  mm = melt(tempS.bwAni)

  new.temp = mm[which(!is.na(match(mm$variable, c("REML.EDF", "LC.EDF", "REAL.EDF")))),]

  new.temp$value = round(new.temp$value, digits = 4)
  colnames(new.temp)[4] = "EDF"


   suppressWarnings(trellis.par.set(layout.heights = list(strip = 1.25),
          superpose.line = list(col=c("black", "red", "blue"), lwd=2,lty=1:3)))     
                                                                                      
   interval = as.numeric(levels(new.temp$gamma.ani)[seq(1,nlevels(new.temp$gamma.ani),2)])
   
   new.temp$gamma.ani = as.numeric(as.character(new.temp$gamma.ani))
   
   xyplot(value ~ gamma.ani|gamma.cag, group = EDF, data = new.temp,lty=1:3, type="l",
   col=c("black", "red", "blue"), lwd=2,
   scales=list(x=list(at=interval, labels=interval, log=TRUE)),
    xlab=expression(paste(sigma[A]^2," / ",sigma^2)),
      key  =
      simpleKey(text = paste(c("Method:  ","", ""), levels(new.temp$EDF)),
          points = FALSE, lines = TRUE, cex = 1, columns=3),
     strip = function(factor.levels, par.strip.text, ...)
    {strip.default(factor.levels =  
    sapply(levels(new.temp$gamma.c), 
    function(x) as.expression(substitute(list(sigma[B]^2/sigma^2==x), list(x=x) ))),
  par.strip.text = trellis.par.get("layout.heights"), ...)},   ... )
   

}
                   
                   

plotEDF2 =
function(sim, rowNames = "Within RunWithinCag.AniResidual", ...){
  require(reshape)
  require(lattice)
  
  tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)),]

  mm = melt(tempS.bwAni)

  new.temp = mm[which(!is.na(match(mm$variable, c("REML.EDF", "LC.EDF", "REAL.EDF")))),]

  new.temp$value = round(new.temp$value, digits = 4)
  colnames(new.temp)[4] = "EDF"


   suppressWarnings(trellis.par.set(layout.heights = list(strip = 1.25),
          superpose.line = list(col=c("black", "red", "blue"), lwd=2,lty=1:3)))     
                                                                                      
   interval = as.numeric(levels(new.temp$gamma.ani)[seq(1,nlevels(new.temp$gamma.ani),2)])
   
   new.temp$gamma.ani = as.numeric(as.character(new.temp$gamma.ani))
   
   xyplot(value ~ gamma.ani|gamma.run, group = EDF, data = new.temp,lty=1:3, type="l",
   col=c("black", "red", "blue"), lwd=2,
   scales=list(x=list(at=interval, labels=interval, log=TRUE)),
    xlab=expression(paste(sigma[A]^2," / ",sigma^2)),
      key  =
      simpleKey(text = paste(c("Method:  ","", ""), levels(new.temp$EDF)),
          points = FALSE, lines = TRUE, cex = 1, columns=3),
     strip = function(factor.levels, par.strip.text, ...)
    {strip.default(factor.levels =  
    sapply(levels(new.temp$gamma.r), 
    function(x) as.expression(substitute(list(sigma[R]^2/sigma^2==x), list(x=x) ))),
  par.strip.text = trellis.par.get("layout.heights"), ...)},   ... )
   

}
                             
                   
                   
                   
                   
                   
                   
                   