incltxt <- '
  arma::umat swap(int s1, int s2, arma::umat X) {

    arma::umat newX = X;
    arma::umat newZ = arma::umat(1,1).zeros();

    arma::urowvec X1 = X.row(s1);
    arma::urowvec X2 = X.row(s2);
            
    if(sum(X1 == X2) == X.n_cols) return(newZ);

    newX.row(s1) = X.row(s2);
    newX.row(s2) = X.row(s1);

    return(newX);
    }
    
  arma::umat check(double aveEff, int oldV, double tr, double trSq, int oldLen,
          double rep, arma::umat newX, arma::mat Y, arma::mat C) {

    arma::umat newZ = arma::umat(1,1).zeros();
    arma::umat newY = arma::umat(1,1).ones();

    arma::mat infoMat = arma::mat(C * arma::trans(newX) * Y * newX * C);
    arma::vec eigval = arma::eig_sym(infoMat);
    
    if(sum(eigval) < 0.0000000001) return(newZ);
      
    arma::vec conEff = eigval.elem(find(eigval > 0.000001))/rep;

    int len = conEff.n_elem;

    //if(len < oldLen) return(newZ);

    double newTr =  arma::trace(infoMat);
    /*
    if( newTr > tr){
      return(newX);

    } else if (newTr == tr){
      */
      double compareAveEff = 1/mean(1/conEff);

      if(compareAveEff > aveEff){
        return(newX);
      }else if(compareAveEff == aveEff){
        arma::mat squInfoMat = square(infoMat);
        double newTrSq = arma::trace(squInfoMat);

        if(newTrSq < trSq){
          return(newX);
        } else if (newTrSq == trSq){

          double newV = var(conEff);

          if( newV < oldV){
            return(newX);
         } else if (newV == oldV){
            return(newY);
         } else {
            return(newZ);
          }
        } else {
          return(newZ);
        }
      } else {
          return(newZ);
      }
    /*
    }  else {
          return(newZ);
    } */
    }
    
    arma::umat lastCheck(double aveEff, int oldV, double tr, double trSq, int oldLen,
          double rep, arma::umat newX, arma::mat Y, arma::mat C) {
        
        
    arma::umat newZ = arma::umat(1,1).zeros();
 
    arma::mat infoMat = arma::mat(C * arma::trans(newX) * Y * newX * C);
        
    arma::vec eigval = arma::eig_sym(infoMat);
    
    if(sum(eigval) < 0.0000000001) return(newZ);
  

    arma::vec conEff = eigval.elem(find(eigval > 0.000001))/rep;

    int len = conEff.n_elem;
      
    //if(len < oldLen) return(newZ);

    double newTr =  arma::trace(infoMat);
    /*
    if( newTr > tr){
      return(newX);

    } else if (newTr == tr){
      */
      double compareAveEff = 1/mean(1/conEff);

      if(compareAveEff > aveEff){
                 
        return(newX);
      }else if(compareAveEff == aveEff){
        arma::mat squInfoMat = square(infoMat);
        double newTrSq = arma::trace(squInfoMat);

        if(newTrSq < trSq){
          return(newX);
        } else if (newTrSq == trSq){

          double newV = var(conEff);

          if( newV < oldV){
            return(newX);
         } else {
            return(newZ);
          }
        } else {
          return(newZ);
        }
      } else {
          return(newZ);
      }
     /*
    }  else {
          return(newZ);
    }
       */
    }
'

 
##############################################################
  

code<- '
  // copy the data to armadillo structures
  arma::umat X = Rcpp::as<arma::umat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat C = Rcpp::as<arma::mat>(C_);
  double rep = as<double>(rep_);
  int tol = as<int>(tol_);

  int nCol = C.n_cols;
  int nRow = C.n_rows;

  int nC = nRow/nCol;

  arma::cube cMat(nCol, nCol, nC) ;
  
  arma::umat newXX, newXX1;
  
  arma::mat infoMat, squInfoMat, eigvec;
  arma::vec eigval, conEff;
  
  arma::vec aveEff(nC), tr(nC), trSq(nC), oldV(nC), oldLen(nC);
  double compareAveEff = 0.0;

  bool exitedInner = false;

  for(int k = 0; k < nC; k++){
   
    cMat.slice(k) = C.rows(k*nCol, (k+1)*nCol -1);
    
    infoMat = arma::mat(cMat.slice(k) * arma::trans(X) * Y * X * cMat.slice(k));
    arma::eig_sym(eigval, eigvec, infoMat, "dc");
    
    if(sum(eigval) < 0.0000000001){
      oldLen(k) = 0;
      aveEff(k) = 0;
      oldV(k) = 1000;
      tr(k) = 0;
      trSq(k) = 0;     
    } else {
      conEff = eigval.elem(find(eigval > 0.000001))/rep;
      
      oldLen(k) = conEff.n_elem;
      aveEff(k) = 1/(mean(1/conEff));
      oldV(k) = var(conEff);
      tr(k) = arma::trace(infoMat);
      squInfoMat = square(infoMat);
      trSq(k) = arma::trace(squInfoMat); 
    }   
  }
 
       
  arma::umat Swaps = arma::umat(1, 2);
  arma::umat swaps = arma::umat(1, 2);

  Swaps.zeros();
  swaps.zeros();
                               
  int xLength = X.n_rows;
  int k;
  int counter = 0;   
    
  for(k = 0; k< tol; k++){
    compareAveEff = aveEff(0);
   
    for( int i = 0; i < (xLength - 1); i++){

      for(int j = i + 1; j < xLength ; j++){
        newXX = swap(i, j, X);  //perform swapping
      
      
        if(newXX.n_rows == 1) continue;

        //First check on C
        if(nC > 1){
          newXX1 = check(aveEff(0), oldV(0), tr(0), trSq(0), oldLen(0), 
                                    rep, newXX, Y, cMat.slice(0));
        } else {
      
          newXX = lastCheck(aveEff(0), oldV(0), tr(0), trSq(0), oldLen(0), 
                                    rep, newXX, Y, cMat.slice(0));
        }         
             
        if((newXX1.n_rows == 1) && (newXX1(0,0) == 0)) continue;
        
        if(newXX.n_rows == 1) continue;
        
   
        for(int l = 1; l < nC; l++){       
          if(l != (nC - 1)){
            if((newXX1.n_rows == 1) && (newXX1(0,0) == 0)){
              exitedInner = true;
              break;
            } else if((newXX1.n_rows == 1) && (newXX1(0,0) == 1)) {
              newXX1 = check(aveEff(l), oldV(l), tr(l), trSq(l), oldLen(l), 
                                  rep, newXX, Y, cMat.slice(l));
            }
          } else {
            //Final check 
            if((newXX1.n_rows == 1) && (newXX1(0,0) == 0)){
              exitedInner = true;
              break;
            } else if((newXX1.n_rows == 1) && (newXX1(0,0) == 1)) {
              newXX = lastCheck(aveEff(l), oldV(l), tr(l), trSq(l), oldLen(l), 
                                  rep, newXX, Y, cMat.slice(l));
            }
          }
        }
         
        if(exitedInner){
          exitedInner = false;
          continue;
        }
        
        if(newXX.n_rows != 1){
        
          swaps(0,0) = i + 1;
          swaps(0,1) = j + 1;
          Swaps.insert_rows(counter, swaps);
          
          counter++;
          X = newXX;
          
          for(int m = 0; m < nC; m++){
            infoMat = arma::mat(cMat.slice(m) * arma::trans(X) * Y * X * cMat.slice(m));
          
            arma::eig_sym(eigval, eigvec, infoMat, "dc");
         
         if(sum(eigval) < 0.0000000001){
              oldLen(m) = 0;
              aveEff(m) = 0;
              oldV(m) = 1000;
              tr(m) = 0;
              trSq(m) = 0;     
            } else {
              conEff = eigval.elem(find(eigval > 0.000001))/rep;
              
              oldLen(m) = conEff.n_elem;
              aveEff(m) = 1/(mean(1/conEff));
              oldV(m) = var(conEff);
              tr(m) = arma::trace(infoMat);
              squInfoMat = square(infoMat);
              trSq(m) = arma::trace(squInfoMat); 
            }   
    
          }
        }
      } 
    }
    if(compareAveEff == aveEff(0)) break;
  }
  
  return Rcpp::List::create(
        Rcpp::Named("X.trt") = X,
        Rcpp::Named("cycle") = k,
        Rcpp::Named("nCan")= oldLen,
        Rcpp::Named("maxAveEff")= aveEff,
        Rcpp::Named("vCan")= oldV,
        Rcpp::Named("swap")= Swaps);
  '


cpp.new.fact.s.opt <- cxxfunction(
 signature(X_ = "matrix", Y_ = "matrix", C_ = "matrix",
 rep_ = "double", tol_ = "int"),
  body = code,  incl = incltxt,plugin = "RcppArmadillo")
 
################################################################################
cpp.new.fact.s.optimised = function(X.trt, blk.proj, C.trt.mat = diag(ncol(X.trt)), 
    tol = 20, nIter = 100) {
    cat("Finding the optimal design using C++ swapping method:\n")
    
    # orig.X.array = array(0, dim = c(nrow(X.trt), ncol(X.trt), nIter)) X.array =
    # array(0, dim = c(nrow(X.trt), ncol(X.trt), nIter)) Swap = vector(mode =
    # 'list', length = nIter)
    
    best.X.trt = X.trt
    orig.X.trt = X.trt
    nC = nrow(C.trt.mat)/ncol(C.trt.mat)
    Swap = matrix()
    
    nCan = aveEff = numeric(nC)
    vCan = rep(Inf, nC)
    
    Rep = nrow(X.trt)/ncol(X.trt)
    
    pb <- txtProgressBar(min = 0, max = nIter, style = 3)
    for (i in 1:nIter) {
        setTxtProgressBar(pb, i)
        
        old.X.trt = X.trt = X.trt[sample(nrow(X.trt)), ]
        XX = cpp.new.fact.s.opt(X.trt, blk.proj, C.trt.mat, Rep, tol)
        # X.trt= X.trt[sample(nrow(X.trt)),]
        
        # XX orig.X.array[,,i] = X.trt X.array[,,i] = XX$X.trt Swap[[i]] =
        # XX$swap[-nrow(XX$swap),]
        
        for (j in 1:nC) {
            
            if (XX$maxAveEff[j, 1] > aveEff[j]) {
                orig.X.trt = old.X.trt
                best.X.trt = XX$X.trt
                nCan = as.numeric(XX$nCan)
                aveEff = as.numeric(XX$maxAveEff)
                vCan = as.numeric(XX$vCan)
                Swap = XX$swap[-nrow(XX$swap), ]
            } else if (XX$maxAveEff[j, 1] == aveEff[j]) {
                if (XX$vCan[j, 1] < vCan[j]) {
                  orig.X.trt = old.X.trt
                  best.X.trt = XX$X.trt
                  nCan = as.numeric(XX$nCan)
                  aveEff = as.numeric(XX$maxAveEff)
                  vCan = as.numeric(XX$vCan)
                  Swap = XX$swap[-nrow(XX$swap), ]
                }
            }
            
        }
    }
    close(pb)
    
    return(list(orig.X.trt = orig.X.trt, best.X.trt = best.X.trt, nCan = nCan, aveEff = aveEff, 
        vCan = vCan, Swap = Swap))
} 

      if (XX$nCan[j, 1] > nCan[j]) {
                orig.X.trt = old.X.trt
                best.X.trt = XX$X.trt
                nCan = as.numeric(XX$nCan)
                aveEff = as.numeric(XX$maxAveEff)
                vCan = as.numeric(XX$vCan)
                Swap = XX$swap[-nrow(XX$swap), ]
                
            } else if (XX$nCan[j, 1] == nCan[1]) {
             }

reSwap = function(des.mat, swap) {
    
    for (i in 1:nrow(swap)) {
        
        temp1 = des.mat[swap[i, 1], ]
        temp2 = des.mat[swap[i, 2], ]
        
        des.mat[swap[i, 2], ] = temp1
        des.mat[swap[i, 1], ] = temp2
        
    }
    
    return(des.mat)
}

 






