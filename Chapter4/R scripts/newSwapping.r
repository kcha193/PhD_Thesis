



swap.stage1.new <- 
function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {

    idx <- seq(1, length(ind) - 1)
    changepoint1 <- sample(idx, size = 1)

    maxC = ceiling(changepoint1/nPlot) * nPlot
    s = unique((maxC-nPlot + 1):maxC)
    
    if(length(s[-which(s == changepoint1)]) == 1){
       changepoint2 = s[-which(s == changepoint1)]
    } else {
      changepoint2 = sample(s[-which(s == changepoint1)], size=1)
    }
      
    changepoint1 = which(Z1.mat[,which(Z1.mat[ind[changepoint1],] == 1)]==1)  
      
    changepoint2 = which(Z1.mat[,which(Z1.mat[ind[changepoint2],] == 1)]==1)  
                        
    tmp <- ind[which(ind == changepoint1)]
    ind[which(ind == changepoint1)] <- ind[which(ind == changepoint2)]
    ind[which(ind == changepoint2)] <- tmp
    return(ind)
}


swap.stage2.new <- 
function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {

    idx <- seq(1, length(ind) - 1)
    changepoint1 <- sample(idx, size = 1)
   
    s = unique(c(seq(changepoint1, max(idx), nPlot), seq(changepoint1, 1, -nPlot)))

    if(length(s[-which(s == changepoint1)]) == 1){
      changepoint2 = s[-which(s == changepoint1)]
    } else {
      changepoint2 = sample(s[-which(s == changepoint1)], size=1)
    }
    
    changepoint1 = which(Z1.mat[,which(Z1.mat[ind[changepoint1],] == 1)]==1)  
      
    changepoint2 = which(Z1.mat[,which(Z1.mat[ind[changepoint2],] == 1)]==1)  
                        
    tmp <- ind[which(ind == changepoint1)]
    ind[which(ind == changepoint1)] <- ind[which(ind == changepoint2)]
    ind[which(ind == changepoint2)] <- tmp
    return(ind)
}

swap.stage3.new <- function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {

    idx <- seq(1, length(ind) - 1)
    changepoint1 <- sample(idx, size = 1)
    #Two checks for omitting the swap of two indtical animals and animals in the same block
    pair.changepoint1 <- getPair(changepoint1, length(ind) - 1)
    
    pair.changepoint1 <- pair.changepoint1[-which(pair.changepoint1==changepoint1)]
    
    maxC = ceiling(changepoint1/nPlot) * nPlot
    pair.maxC = ceiling(pair.changepoint1/nPlot) * nPlot

    if(pair.changepoint1 == 1){
      toInitial =  pair.changepoint1
    } else{
      toInitial = seq(pair.changepoint1, 1, -nPlot)
    }

    if(pair.changepoint1 == (length(ind) -1)){
      toFinal=  pair.changepoint1
    } else{
      toFinal = seq(pair.changepoint1, length(ind) -1, nPlot)
    }

    s = unique(c((maxC-nPlot + 1):maxC, seq(changepoint1, length(ind) -1 , nPlot), 
                                        seq(changepoint1, 1, -nPlot),
                 (pair.maxC-nPlot + 1):pair.maxC,  toInitial, toFinal))
    
    if(length((1:(length(ind)-1))[-s])<1){
      return(ind)
    } else{
      changepoint2 = sample((1:(length(ind)-1))[-s], size = 1)
    }
    
    tmp <- ind[getPair(changepoint1, length(ind) )]
    ind[getPair(changepoint1, length(ind) )] <- ind[getPair(changepoint2, length(ind))]
    ind[getPair(changepoint2, length(ind))] <- tmp
    return(ind)
}
