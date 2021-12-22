

jumpState =
function(X.trt, blk.proj){



  states = lapply(1:nrow(X.trt), function(x) t(cbind(x, (x+1):nrow(X.trt))))
  states = states[-length(states)]

  #combine all the states into a matrix
  states= t(matrix(c(states, recursive=TRUE), nrow = 2))

  swap = apply(states, 1, function(x) update.Infor.mat(x, X.trt, blk.proj))

  if(is.null(swap)) return(X.trt)
  
  Aa1 = c(lapply(swap, function(x) ifelse(is.null(x$Aa), NA, x$Aa)), recursive=TRUE)

  bestSwap1 = which(Aa1 == max(Aa1, na.rm = TRUE))

  Ab1 = swap[[bestSwap1[1]]]$AbAveEff

  if(Ab1 > max(Aa1, na.rm = TRUE)) return(X.trt)

  X.trt1 = lapply(swap[bestSwap1], function(x) x$newTrtMat)

  for(l in 1:nrow(states)){
    #print(l)
    swap2 = apply(states, 1, function(x) lapply(X.trt1, function(y) update.Infor.mat(x, y, blk.proj)))

    #Break the list of list to list
    swap2.list = list()
    for( i in 1:length(swap2)){
      s = swap2[[i]]
      for(j in 1:length(s)){
        if(is.null(s[[j]])) next

        swap2.list = c(swap2.list, s[j])
      }
    }

    swap2 = swap2.list

    if(length(swap2) ==0) return(X.trt1)

    Aa2 = c(lapply(swap2, function(x) ifelse(is.null(x$Aa), NA, x$Aa)), recursive=TRUE)

    bestSwap2 = which(Aa2 == max(Aa2, na.rm = TRUE))

    X.trt1 = lapply(swap2[bestSwap2], function(x) x$newTrtMat)
    
    index = sapply(X.trt1, function(x) sapply( X.trt1, function(y) all(y==x)))

    index[lower.tri(index, diag = TRUE)] = FALSE
    
    X.trt1 = X.trt1[!apply(index, 2, any)]
    

    if(length(X.trt1) > 20){

      index = sample(20)
      X.trt1 = X.trt1[as.numeric(index[1:(20*0.5)])]

    }

  }

}

