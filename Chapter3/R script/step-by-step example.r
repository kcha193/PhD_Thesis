

extractAniTrt = function(bestInd){
  new.Z1.mat=  Z1.mat[bestInd,]
  #new.Z1.mat=  Z1.mat

  ################################################################################
  #Set the design from the search

  colnames(new.Z1.mat) = sort(levels(interaction(multiLetters(1:nZ1))))
  new.Z1.des = apply(new.Z1.mat, 1,  function(x)
  colnames(new.Z1.mat)[which(as.logical(x))])
  new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))

  design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

design.df = cbind(design.df,
            Ani = t(new.Z1.des),
            Trt = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                          as.character(phase1DesignEX1$Ani)),]$Trt)

  print(matrix(design.df$Ani, nrow = nBlk, ncol = nPlot, byrow = TRUE))
  print(matrix(tolower(design.df$Trt), nrow = nBlk, ncol = nPlot, byrow = TRUE))

}


print(c(nTrt, bRep, tRep, nPlot))
   
newInit = c(1:n, 1)
extractAniTrt(1:n)

(old = obj.fun.CRD(c(1:n, 1)) )

iter = 3e4

cat("Level: 1, Finding the temperatures for both stages\n")

temp = initialTemp(iter, newInit, swap.stage2.new, obj.fun.CRD)

res.stage2 = temp$newInit
y2.stage2 = temp$y2
y2low.stage2 = temp$y2low


temp = initialTemp(iter, newInit, swap.stage1.new, obj.fun.CRD)

res.stage1 = temp$newInit
y2.stage1 = temp$y2
y2low.stage1 = temp$y2low


temp = initialTemp(iter, newInit,swap.stage3.new, obj.fun.CRD)

res = temp$newInit
y2 = temp$y2
y2low = temp$y2low

cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
cat("Stage 0 temperature range: ", y2, " to ", y2low, "\n")

obj.fun.new = obj.fun.CRD 
 obj.fun = obj.fun.CRD 
 test.obj.fun.new = test.obj.fun.CRD 
 
test.newInit = apply(rbind(newInit, res,res.stage1,res.stage2),1, obj.fun.new)
newInit = rbind(newInit, res,res.stage1,res.stage2)[which(test.newInit == max(test.newInit))[1],]

(cur = obj.fun.new(newInit))

extractAniTrt(newInit[-length(newInit)])

while ((cur - old) > 1e-06) {
      
      temp = initialTemp(iter, newInit, swap.stage2.new, obj.fun)
  
      res.stage2 = temp$newInit
      y2.stage2 = temp$y2
      y2low.stage2 = temp$y2low
  
  
      temp = initialTemp(iter, newInit, swap.stage1.new, obj.fun)
  
      res.stage1 = temp$newInit
      y2.stage1 = temp$y2
      y2low.stage1 = temp$y2low
  
  
      temp = initialTemp(iter, newInit, swap.stage3.new, obj.fun)
  
      res = temp$newInit
      y2 = temp$y2
      y2low = temp$y2low
         
      test.newInit = apply(rbind(res,res.stage1,res.stage2),1, obj.fun)
      newInit = rbind(res,res.stage1,res.stage2)[which(test.newInit == max(test.newInit))[1],]

      old = cur
      
      cur = obj.fun(newInit)
      print(test.obj.fun.new(newInit))     
 
      cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
      cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
      cat("Stage 0 temperature range: ", y2, " to ", y2low, "\n")
  
}


bestInit =  newInit
bestObj = obj.fun.new(newInit)
mem.bestInit = numeric(3e4)
mem.newInit = numeric(3e4)

first.newInit = newInit

y = y2.stage2 / log((((1:iter)-1) %/% 1000)*1000 + exp(1))


################################################################################
extractAniTrt(1:n)

extractAniTrt(newInit[-length(newInit)])

(x1 = obj.fun.new(newInit) )

#First swap of the first stage
newInit1 = swap.stage1.new(newInit)
extractAniTrt(newInit1[-length(newInit1)])


obj.fun.new(newInit1)

#Check for based on the simulated annealing Temperture of 100

exp((obj.fun.new(newInit1)  - x1)/y[1])

set.seed(1650)
(acceptancy.prop = runif(1,0,1))

if(exp((obj.fun.new(newInit1)  - x1)/y[1]) >  acceptancy.prop){
  print("improved 1!")
  mem.newInit[1] = obj.fun.new(newInit1)
}

if((obj.fun.new(newInit1) > bestObj)){
  print("improved 2!")

  bestObj =  obj.fun.new(newInit1)
  bestInit =  newInit1
  mem.bestInit[1] = bestObj
          print(bestObj)

}

#Second swap
newInit2 = swap.stage1.new(newInit1)
extractAniTrt(newInit2[-length(newInit2)])


obj.fun.new(newInit2)

#Check for based on the simulated annealing Temperture of 100

exp((obj.fun.new(newInit2)  - obj.fun.new(newInit1))/y[2])

set.seed(1704)
(acceptancy.prop = runif(1,0,1))

if(exp((obj.fun.new(newInit2)  - obj.fun.new(newInit1))/y[2]) >  acceptancy.prop){
  print("improved 1!")
    mem.newInit[2] = obj.fun.new(newInit2)
}

if((obj.fun.new(newInit2) > bestObj)){
 print("improved 2!")

  bestObj =  obj.fun.new(newInit2)
  bestInit =  newInit2
    mem.bestInit[2] = bestObj
          print(bestObj)

}

#Third swap
newInit3 = swap.stage1.new(newInit2)
extractAniTrt(newInit3[-length(newInit3)])


obj.fun.new(newInit3)

#Check for based on the simulated annealing Temperture of 100

exp((obj.fun.new(newInit3)  - obj.fun.new(newInit2))/y[3])

(acceptancy.prop = runif(1,0,1))

if(exp((obj.fun.new(newInit3)  - obj.fun.new(newInit2))/y[3]) >  acceptancy.prop){
  print("improved 1!")
      mem.newInit[3] = obj.fun.new(newInit3)

}

if((obj.fun.new(newInit3) > bestObj)){
 print("improved 2!")

  bestObj =  obj.fun.new(newInit3)
  bestInit =  newInit3
      mem.bestInit[3] = bestObj
          print(bestObj)

}

#Third swap

for(i in 4:(3e4)){
  newInit = swap.stage1.new(newInit3)

  #cat("Step", i, " \n")

  #print(obj.fun.new(newInit))
  #extractAniTrt(newInit[-length(newInit)])

  (acceptancy.prop = runif(1,0,1))

  #Check for based on the simulated annealing Temperture of 100

  #print(exp((obj.fun.new(newInit)  - obj.fun.new(newInit3))/y[i]))


  if(exp((obj.fun.new(newInit)  - obj.fun.new(newInit3))/y[i]) >  acceptancy.prop){
    newInit3 = newInit
    #print("improved 1!")
          mem.newInit[i] = obj.fun.new(newInit)
  }

  if((obj.fun.new(newInit) > bestObj)){
   print("improved 2!")

    bestObj =  obj.fun.new(newInit)
    bestInit =  newInit
    #print(bestObj)
          mem.bestInit[i] = bestObj

  }
  
  
}

bestObj
extractAniTrt(bestInit[-length(bestInit)])

obj.fun.new(c(1:n,1))

obj.fun.new(bestInit)

mem.bestInit

plot(1:3e4, mem.newInit[1:3e4], type = "l")

plot(1:3e4, mem.bestInit[1:3e4], type = "l")


####################################################################################################
#First swap of the second stage

y = y2.stage1 / log((((1:30000)-1) %/% 1000)*1000 + exp(1))


newInit3 = swap.stage2.new(bestInit)

 print(obj.fun.new(newInit3))
  extractAniTrt(newInit3[-length(newInit3)])

  #Check for based on the simulated annealing Temperture of 100
   set.seed(1928)
   (acceptancy.prop = runif(1,0,1))

  print(exp((obj.fun.new(newInit3)  - obj.fun.new(bestInit))/y[1]))


  if(exp((obj.fun.new(newInit3)  - obj.fun.new(bestInit))/y[1]) >  acceptancy.prop){
    newInit3 = newInit3
    print("improved 1!")
  }

  if((obj.fun.new(newInit3) > bestObj)){
   print("improved 2!")

    bestObj =  obj.fun.new(newInit3)
    bestInit =  newInit
    print(bestObj)

  }


#Third swap
for(i in 2:(3e4)){
  newInit = swap.stage2.new(newInit3)

  #cat("Step", i, " \n")

  #print(obj.fun.new(newInit))
  #extractAniTrt(newInit[-length(newInit)])

  #Check for based on the simulated annealing Temperture of 100
   (acceptancy.prop = runif(1,0,1))

  #print(exp((obj.fun.new(newInit)  - obj.fun.new(newInit3))/y[i]))


  if(exp((obj.fun.new(newInit)  - obj.fun.new(newInit3))/y[i]) >  acceptancy.prop){
    newInit3 = newInit
    #print("improved 1!")
  }

  if((obj.fun.new(newInit) > bestObj)){
   print("improved 2!")

    bestObj =  obj.fun.new(newInit)
    bestInit =  newInit
     #   print(bestObj)

  }
}


bestObj
extractAniTrt(bestInit[-length(bestInit)])

obj.fun.new(bestInit)

################################################################################

####################################################################################################
#First swap of the third stage
y = y2 / log((((1:30000)-1) %/% 1000)*1000 + exp(1))

newInit3 = swap.stage3.new(bestInit)

 print(obj.fun.new(newInit3))
  extractAniTrt(newInit3[-length(newInit3)])

  print(exp((obj.fun.new(newInit3)  - obj.fun.new(bestInit))/y[1]))
  
  #Check for based on the simulated annealing Temperture of 100
  set.seed(44)
   (acceptancy.prop = runif(1,0,1))

  #print(exp((obj.fun.new(newInit3)  - obj.fun.new(bestInit))/y))


  if(exp((obj.fun.new(newInit3)  - obj.fun.new(bestInit))/y[1]) >  acceptancy.prop){
    newInit3 = newInit3
    print("improved 1!")
  }

  if((obj.fun.new(newInit3) > bestObj)){
   print("improved 2!")

    bestObj =  obj.fun.new(newInit3)
    bestInit =  newInit
    print(bestObj)

  }


#Third swap
for(i in 2:(3e4)){
  newInit = swap.stage3.new(newInit3)

  #cat("Step", i, " \n")

 # print(obj.fun.new(newInit))
 # extractAniTrt(newInit[-length(newInit)])

  #Check for based on the simulated annealing Temperture of 100
   (acceptancy.prop = runif(1,0,1))

 # print(exp((obj.fun.new(newInit)  - obj.fun.new(newInit3))/y[i]))


  if(exp((obj.fun.new(newInit)  - obj.fun.new(newInit3))/y[i]) >  acceptancy.prop){
    newInit3 = newInit
  #  print("improved 1!")
  }

  if((obj.fun.new(newInit) > bestObj)){
  # print("improved 2!")

    bestObj =  obj.fun.new(newInit)
    bestInit =  newInit
        print(bestObj)
 #
  }
}


bestObj
extractAniTrt(bestInit[-length(bestInit)])

obj.fun.new(bestInit)

################################################################################
