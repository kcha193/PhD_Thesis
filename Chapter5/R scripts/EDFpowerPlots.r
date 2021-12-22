


 colnames(temp1) = c("gamma.run", "gamma.ani", "DF", "Initial Power", "EDF", "Adjusted Power")
 
  temp1 = remainNeg

  temp1$gamma.run = as.factor(temp1$gamma.run)
   temp1$gamma.ani = as.factor(temp1$gamma.ani)

  temp1$diffNeg =  as.numeric(as.character(temp1$"Adjusted MS")) - 
                as.numeric(as.character(temp1$MS))
  
 
        mm = melt(temp1[,c("gamma.run", "gamma.ani", "Initial Power","Adjusted Power")])

    #names(mm)[4] = "Method"

        names(mm)[3] = "Adjustment"


    mm$gamma.ani = as.numeric(as.character(mm$gamma.ani))
    # mm$gamma.run = as.factor( paste('gamma.run =', mm$gamma.run)) levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))


    pg <- ggplot(mm, aes(x = log(gamma.ani), y = value)) +
    ylab("F-test Power") + xlab(expression(paste(sigma[A]^2, " / ", sigma^2)))  +  ggtitle("Power plot") + 
    scale_x_continuous(breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)],
            labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) +
    geom_line(aes(colour = Adjustment ), size = 1, alpha = I(0.5)) +
    facet_wrap(~gamma.run)  + theme(legend.position="top")

   pg


  pdf("remainNegPowerEDF.pdf", width = 12)
    pg
  dev.off()

#####################################################################################################################    
   
  temp1 = adjustZero

#   colnames(temp1) = c("gamma.run", "gamma.ani", "DF", "Initial Power", "EDF", "Adjusted Power")
 
   temp1$gamma.run = as.factor(temp1$gamma.run)
   temp1$gamma.ani = as.factor(temp1$gamma.ani)

  
  mm = melt(temp1[,c("gamma.run", "gamma.ani", "Initial Power", "Adjusted Power")])

 
    names(mm)[3] = "Adjustment"

    mm$gamma.ani = as.numeric(as.character(mm$gamma.ani))
    # mm$gamma.run = as.factor( paste('gamma.run =', mm$gamma.run)) levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))


      pg <- ggplot(mm, aes(x = log(gamma.ani), y = value)) +
    ylab("F-test Power") + xlab(expression(paste(sigma[A]^2, " / ", sigma^2)))  +  ggtitle("Power plot") + 
    scale_x_continuous(breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)],
            labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) +
    geom_line(aes(colour = Adjustment ), size = 1, alpha = I(0.5)) +
    facet_wrap(~gamma.run)  + theme(legend.position="top")

   pg

  pdf("adjustZeroPowerEDF.pdf", width = 12)
    pg
  dev.off()
   
   
   
   
   
   