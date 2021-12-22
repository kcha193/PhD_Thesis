
incltxt <- '
  arma::umat swap(double aveEff, double oldV, int oldLen,  
          double rep, int s1, int s2, arma::umat X, arma::mat Y, arma::mat C) {
    
    arma::umat newX = X;
    arma::umat newZ = arma::umat(1,1).zeros();    
    
    arma::urowvec X1 = X.row(s1); 
    arma::urowvec X2 = X.row(s2);
    
    if(sum(X1 == X2) == X.n_cols) return(newZ);
    
    newX.row(s1) = X.row(s2);
    newX.row(s2) = X.row(s1);
    
    arma::mat infoMat = arma::mat(C * arma::trans(newX) * Y * newX * C);

    arma::vec eigval = arma::eig_sym(infoMat);
  
    arma::vec conEff = eigval.elem(find(eigval > 0.000001))/rep;
  	
    
    int len = conEff.n_elem;
    
    if(len < oldLen) return(newZ);
    
    double compareAveEff= 1/mean(1/conEff);
     
    if(compareAveEff > aveEff){ 
      return(newX);
    }else if(compareAveEff == aveEff){  
      double newV = var(conEff);    
            
      if( newV< oldV){
        return(newX);      
      } else {
        return (newZ);     
      }     
    }else{
      return (newZ);
    }   
  }
'

code<- '
  // copy the data to armadillo structures
  arma::umat X = Rcpp::as<arma::umat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat Z = Rcpp::as<arma::mat>(Z_);
  arma::mat C = Rcpp::as<arma::mat>(C_);
    double rep = as<double>(rep_);
      int tol = as<int>(tol_);

   double maxAveEff = 0.0;
  arma::mat infoMat;
  arma::vec eigval;
  arma::vec conEff;

  double aveEff;
  int oldLen;

   infoMat = arma::mat(C * arma::trans(X) * Y * X * C);
      eigval = arma::eig_sym(infoMat);
      conEff = eigval.elem(find(eigval > 0.000001))/rep;
      oldLen = conEff.n_elem;
      aveEff = 1/(mean(1/conEff));

        arma::rowvec z = Z.row(tol);
        
        int s1 = z[0];
        int s2 = z[1];
        arma::umat newXX = swap(aveEff, oldLen, rep, s1, s2, X, Y, C);
        
        if(newXX.n_rows != 1){
            
          infoMat = arma::mat(C * arma::trans(newXX) * Y * newXX * C);
          eigval = arma::eig_sym(infoMat);
          conEff = eigval.elem(find(eigval > 0.000001))/rep;
               
          oldLen = conEff.n_elem;               
          aveEff = 1/mean(1/conEff);
          
          if(maxAveEff < aveEff) maxAveEff = aveEff;
        }
           
    return Rcpp::wrap(aveEff);
  '
  
cpp.test <- cxxfunction(
 signature(X_ = "matrix", Y_ = "matrix", Z_ = "matrix", C_ = "matrix", 
 rep_ = "double", tol_ = "int"),
  body = code, incl = incltxt, plugin = "RcppArmadillo")
  
cpp.test(X.trt, blk.proj, cpp.states, C.trt.mat, Rep,3)   
 
  
code<- '
  // copy the data to armadillo structures
  arma::umat X = Rcpp::as<arma::umat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat Z = Rcpp::as<arma::mat>(Z_);
  arma::mat C = Rcpp::as<arma::mat>(C_);
  double rep = as<double>(rep_);
  int tol = as<int>(tol_);
  arma::ucube xMats(X.n_rows, X.n_cols, tol);
  arma::ucube xMats1(X.n_rows, X.n_cols, tol);  
  arma::ucube xMats2(X.n_rows, X.n_cols, tol);

  double maxAveEff = 0.0;
  arma::mat infoMat;
  arma::vec eigval;
  arma::vec conEff;

  double aveEff;
  double oldV;
  
  double newAveEff;
  int oldLen;
  
  xMats.zeros();
  xMats1.zeros();
  
  arma::umat newXX;
    
  int length = Z.n_rows; 
  
  xMats.slice(0) = X; 
  
  int lenList = 1;
  
  int counter = 0;

  for( int k = 0; k < length; k++){
    
    for( int j = 0; j < lenList; j++){
      if(counter >= tol) break;  
      
      X = xMats.slice(j);
      
      infoMat = arma::mat(C * arma::trans(X) * Y * X * C);
      eigval = arma::eig_sym(infoMat);
      conEff = eigval.elem(find(eigval > 0.000001))/rep;
      oldLen = conEff.n_elem;
      aveEff = 1/(mean(1/conEff));
      oldV = var(conEff);
        
      for( int i = 0; i < length; i++){
        if(counter >= tol) break;
                
        arma::rowvec z = Z.row(i);
        
        int s1 = z[0];
        int s2 = z[1];
        
        newXX = swap(aveEff, oldV, oldLen, rep, s1, s2, X, Y, C);
         
        //get the globel max  average efficiency factor      
        if(newXX.n_rows != 1){
          xMats1.slice(counter) = newXX;
          counter++;
          
          infoMat = arma::mat(C * arma::trans(newXX) * Y * newXX * C);
          eigval = arma::eig_sym(infoMat);
          conEff = eigval.elem(find(eigval > 0.000001))/rep;
               
          oldLen = conEff.n_elem;               
          newAveEff = 1/mean(1/conEff);
          
          if(maxAveEff < newAveEff) maxAveEff = newAveEff;
        }
      }
    }       
    
    if(counter == 0) break;         
    
    xMats2.zeros();
    lenList = counter;
    
    counter = 0;
    
    for( int l = 0; l < lenList; l++){
      X = xMats1.slice(l);
      
      infoMat = arma::mat(C * arma::trans(X) * Y * X * C);
      eigval = arma::eig_sym(infoMat);
      conEff = eigval.elem(find(eigval > 0.000001) )/rep;
     
      int len = conEff.n_elem;
      
      if(len < oldLen) continue;
        
      aveEff =  1/mean(1/conEff);
  
      if(maxAveEff == aveEff){
        xMats2.slice(counter) = X;
        counter++;
      }
    }
     
    xMats.zeros();
    xMats = xMats2;
    
    xMats1.zeros();
    
    lenList = counter;
    
    counter = 0;
  }
      
  return Rcpp::List::create( 
  Rcpp::Named("X.trt") = xMats,
    Rcpp::Named("length")= lenList,
     Rcpp::Named("maxAveEff")= maxAveEff);
'


cpp.s.opt <- cxxfunction(
 signature(X_ = "matrix", Y_ = "matrix", Z_ = "matrix", C_ = "matrix", 
 rep_ = "double", tol_ = "int"),
  body = code, incl = incltxt, plugin = "RcppArmadillo")

n.eff = 
function(x, blk.proj, Rep){
  e.va =  (eigen(t(x) %*% blk.proj %*% x)$va)/Rep
  if(all(e.va>1e-7)){
    length(e.va)
  }else {
    length(e.va[-which(e.va<1e-7)])
  }

}


ave.eff=
function(x, blk.proj, Rep){
  e.va =  (eigen(t(x) %*% blk.proj %*% x)$va)/Rep
  if(all(e.va>1e-7)){
     1/mean(1/e.va) 
  } else{
    e.va = e.va[-which(e.va<1e-7)]
    1/mean(1/e.va)
  }
}

bal.eff=
function(x, blk.proj, Rep){
  e.va =  eigen(t(x) %*% blk.proj %*% x)$va/Rep
  if(all(e.va>1e-7)){
    sd(e.va)
  }else {
    sd(e.va[-which(e.va<1e-7)])
  }
}

cpp.s.optimised =
function(X.trt, blk.proj, C.trt.mat = diag(ncol(X.trt)), tol = 100, nIter = 10){
  cat("Finding the optimal design using C++ swapping method:\n")

  X.array = array(0, dim = c(nrow(X.trt), ncol(X.trt), tol*nIter))
  states = lapply(1:nrow(X.trt), function(x) t(cbind(x, (x+1):nrow(X.trt))))
  states = states[-length(states)]
  
  #combine all the states into a matrix
  states= t(matrix(c(states, recursive=TRUE), nrow = 2))
       
  cpp.states = states - 1
  
  Rep = nrow(X.trt)/ncol(X.trt)
  count = 1;
  pb <- txtProgressBar(min = 0, max = nIter, style = 3)
  for( i in 1:nIter){
    setTxtProgressBar(pb, i)
    
    if(nrow(cpp.states) > 500){
      cpp.states = cpp.states[sample(nrow(cpp.states)),]
      cpp.states = cpp.states[1:500,]
    }    
         
    X.trt= X.trt[sample(nrow(X.trt)),]
    XX = cpp.s.opt(X.trt, blk.proj, cpp.states, C.trt.mat, Rep, tol) 
    counter = XX$length
    X.array[,,count:(count+counter-1)] = XX$X.trt[,,1:counter]
    count = count+counter
  }
  close(pb)
  
  X.array = X.array[,,1:(count-1)]
  neff = apply(X.array, 3, function(x) n.eff(x, blk.proj, Rep))
  X.array= X.array[,,which(neff == max(neff, na.rm = TRUE))]

  aveeff = apply(X.array, 3, function(x) ave.eff(x, blk.proj, Rep))
  X.array= X.array[,,which(aveeff == max(aveeff, na.rm = TRUE))]

  if(length(dim(X.array)) == 3){
    baleff = apply(X.array, 3, function(x) bal.eff(x, blk.proj, Rep))
    X.array= X.array[,,which(baleff == min(baleff, na.rm = TRUE))]
  }

  return(X.array)
}

################################################################################
################################################################################
#R functions
jumpState =
function(X.trt, blk.proj, Rep, states){
    
  swap = apply(states, 1, function(x) update.Infor.mat(x, X.trt, blk.proj, Rep))

  if(is.null(swap)) return(X.trt)

  Aa1 = c(lapply(swap, function(x) ifelse(is.null(x$Aa), NA, x$Aa)), recursive=TRUE)

  bestSwap1 = which(Aa1 == max(Aa1, na.rm = TRUE))

  Ab1 = swap[[bestSwap1[1]]]$AbAveEff

  if(Ab1 > max(Aa1, na.rm = TRUE)) return(X.trt)

  X.trt1 = lapply(swap[bestSwap1], function(x) x$newTrtMat)

  for(l in 1:nrow(states)){
    swap2 = apply(states, 1, function(x) lapply(X.trt1, function(y) update.Infor.mat(x, y, blk.proj, Rep)))

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

    index = sapply(X.trt1, function(x) sapply( X.trt1, function(y) all(y == x)))

    index[lower.tri(index, diag = TRUE)] = FALSE

    X.trt1 = X.trt1[!apply(index, 2, any)]


    if(length(X.trt1) > 20){

      index = sample(20)
      X.trt1 = X.trt1[as.numeric(index[1:(20*0.5)])]

    }

  }

}


update.Infor.mat =
function(int, trt.mat, blk.proj, Rep){

  library(MASS)

  Ab = t(trt.mat) %*% blk.proj %*% trt.mat
  ginv.Ab = ginv(Ab)
  Rep = nrow(trt.mat)/ncol(trt.mat)

  a = int[1]
  b = int[2]

  T1 = trt.mat

  x1 = T1[b,] = trt.mat[a,]
  x2 = T1[a,] = trt.mat[b,]

  if(isTRUE(all.equal(x1,x2))) return(NULL)

  c1 = blk.proj[a,a]
  c2 = blk.proj[b,b]

  d1 = blk.proj[a, -c(a,b)]
  d2 = blk.proj[b, -c(a,b)]

  T2 = trt.mat[-c(a,b),]


  B = (x1 - x2) %*% t(d1 - d2) %*% T2

  #return( list( n.trt.mat = T1,
  #A = (t(trt.mat) %*% blk.proj %*% trt.mat) -
  Aa = Ab - (((c1-c2) * (x1%*% t(x1) -x2%*%t(x2))) + B + t(B))


  e.va = eigen(Aa)$va
  e.va = e.va[-which(e.va<1e-7)]
  Aa.ave.eff = 1/mean(Rep/e.va)

  len = length(e.va)

  e.va = eigen(Ab)$va
  e.va = e.va[-which(e.va<1e-7)]
  Ab.ave.eff =  1/mean(Rep/e.va)
  
  if(len < length(e.va)) return(NULL)
  
  if(is.na(Aa.ave.eff)){
    return (NULL)
  }else if(isTRUE(all.equal(Aa.ave.eff, Ab.ave.eff)) ){
     return(NULL) 
  } else if (Aa.ave.eff > Ab.ave.eff){
    return (list(AbAveEff = Ab.ave.eff,
                newTrtMat = T1,
                AaAveEff = Aa.ave.eff))             
  }else{
    return (NULL)
   
  }
}

ave.eff=
function(x, blk.proj, Rep){
  e.va =  eigen(t(x) %*% blk.proj %*% x)$va
  if(all(e.va>1e-7)){
    (length(e.va)/sum(1/e.va))/Rep
  } else{
    (length(e.va[-which(e.va<1e-7)])/sum(1/e.va[-which(e.va<1e-7)]))/Rep
  }
}

bal.eff=
function(x, blk.proj, Rep){
  e.va =  eigen(t(x) %*% blk.proj %*% x)$va

  sd(e.va[-which(e.va<1e-7)])
}



s.optimised=
function(X.trt, blk.proj, nIter = 20){
cat("Finding the optimal design using swapping method:\n")
  
  states = lapply(1:nrow(X.trt), function(x) t(cbind(x, (x+1):nrow(X.trt))))
  states = states[-length(states)]

  #combine all the states into a matrix
  states= t(matrix(c(states, recursive=TRUE), nrow = 2))

  X.trt1 = list()
  Rep = nrow(X.trt)/ncol(X.trt)

  pb <- txtProgressBar(min = 0, max = nIter, style = 3)
  
  for( i in 1:nIter){
    setTxtProgressBar(pb, i)

    X.trt= X.trt[sample(nrow(X.trt)),]  #restart a new design, atempt to jump out of the local max-/minmum
    #swapping two design points, as like moving from state to state
    temp.X.trt = jumpState(X.trt, blk.proj, Rep, states) 
    
    if( i == 1){
      X.trt1 = temp.X.trt
    } else{

      if(is.matrix(temp.X.trt)){
        index = sapply(X.trt1, function(x)  all(temp.X.trt==x))

        X.trt1 = X.trt1[which(!index)]
        X.trt1[[length(X.trt1)+1]] = temp.X.trt
      }else{
        index = sapply(X.trt1, function(x) sapply(temp.X.trt, function(y) all(y==x)))
        if(length(index) == 0){

          X.trt1 = temp.X.trt
        }else{
          X.trt1 = c(X.trt1[which(!apply(index, 2, any))], temp.X.trt)
        }
      }
      ave.eff.X.trt =
        sapply( X.trt1, function(x) ave.eff(x, blk.proj = blk.proj, Rep = Rep))

      X.trt1 =  X.trt1[which(ave.eff.X.trt == max(ave.eff.X.trt))]

      sd.eff.X.trt=
        sapply(X.trt1, function(x) bal.eff(x, blk.proj = blk.proj, Rep = Rep))

      X.trt1 =  X.trt1[which(sd.eff.X.trt == min(sd.eff.X.trt))]

    }
  }
  close(pb)

  return(X.trt1)
}


