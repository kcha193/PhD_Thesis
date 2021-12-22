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

  arma::umat check(double aveEff, double oldV, double tr, double trSq, int oldLen,
          double rep, arma::umat newX, arma::mat Y) {

    arma::umat newZ = arma::umat(1,1).zeros();

    arma::mat infoMat = arma::mat(arma::trans(newX) * Y * newX);
    arma::vec eigval = arma::eig_sym(infoMat);
    arma::vec conEff = eigval.elem(find(eigval > 0.000001))/rep;

    int len = conEff.n_elem;

    if(len < oldLen) return(newZ);

    double newTr =  arma::trace(infoMat);

    if( newTr > tr){
      return(newX);

    } else if (newTr == tr){

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
         }else {
            return(newZ);
          }
        } else {
          return(newZ);
        }
      } else {
          return(newZ);
      }

    }  else {
          return(newZ);
    }

    }
'


code<- '
  // copy the data to armadillo structures
  arma::umat X = Rcpp::as<arma::umat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  double rep = as<double>(rep_);
  int tol = as<int>(tol_);

  arma::umat Swaps = arma::umat(1, 2);
  arma::umat swaps = arma::umat(1, 2);

  Swaps.zeros();
   swaps.zeros();

  arma::mat infoMat, squInfoMat;
  arma::vec eigval;
  arma::vec conEff;

  double aveEff, tr, trSq;
  double compareAveEff = 0.0;
  int oldLen;

  double oldV;

  infoMat = arma::mat(arma::trans(X) * Y * X);
  eigval = arma::eig_sym(infoMat);
  conEff = eigval.elem(find(eigval > 0.000001))/rep;
  oldLen = conEff.n_elem;
  aveEff = 1/(mean(1/conEff));
  oldV = var(conEff);
  tr = arma::trace(infoMat);
  squInfoMat = square(infoMat);
  trSq = arma::trace(squInfoMat);

  int xLength = X.n_rows;
  int k;
  int counter = 0;
  for(k = 0; k< tol; k++){
    //Rprintf("Cycle %d", k);
    //Rprintf("\\n");

    compareAveEff = aveEff;
    for( int i = 0; i < (xLength - 1); i++){

      for(int j = xLength - 1; j >i; j--){

        arma::umat newXX = swap(i, j, X);

        if(newXX.n_rows == 1) continue;

        newXX = check(aveEff, oldV, tr, trSq, oldLen, rep, newXX, Y);

        //std::cout << "value at "<< k << std::endl;
        if(newXX.n_rows != 1){
            swaps(0,0) = i + 1;
            swaps(0,1) = j + 1;
            Swaps.insert_rows(counter, swaps);
            counter++;
            X = newXX;
            infoMat = arma::mat(arma::trans(X) * Y * X );
            eigval = arma::eig_sym(infoMat);
            conEff = eigval.elem(find(eigval > 0.000001))/rep;
            oldLen = conEff.n_elem;
            aveEff = 1/(mean(1/conEff));
            oldV = var(conEff);
            tr = arma::trace(infoMat);
            squInfoMat = square(infoMat);
            trSq = arma::trace(squInfoMat);

        }
      }
    }
    if(compareAveEff == aveEff) break;
  }
  return Rcpp::List::create(
        Rcpp::Named("X.trt") = X,
        Rcpp::Named("cycle") = k,
        Rcpp::Named("maxAveEff")= aveEff,
        Rcpp::Named("swap")= Swaps);
  '


cpp.new.s.opt <- cxxfunction(
 signature(X_ = "matrix", Y_ = "matrix",
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

cpp.new.s.optimised =
function(X.trt, blk.proj, tol = 10, nIter = 20){
  cat("Finding the optimal design using C++ swapping method:\n")

  orig.X.array = array(0, dim = c(nrow(X.trt), ncol(X.trt), nIter))
  X.array = array(0, dim = c(nrow(X.trt), ncol(X.trt), nIter))
  Swap = vector(mode = "list", length = nIter)
  Rep = nrow(X.trt)/ncol(X.trt)

  pb <- txtProgressBar(min = 0, max = nIter, style = 3)
  for( i in 1:nIter){
    setTxtProgressBar(pb, i)

    X.trt= X.trt[sample(nrow(X.trt)),]
    XX = cpp.new.s.opt(X.trt, blk.proj, Rep, tol)

    orig.X.array[,,i] = X.trt
    X.array[,,i] = XX$X.trt
    Swap[[i]] = XX$swap[-nrow(XX$swap),]
  }
  close(pb)

  neff = apply(X.array, 3, function(x) n.eff(x, blk.proj, Rep))
  index = which(neff == max(neff, na.rm = TRUE))
  orig.X.array = orig.X.array[,,index]
  X.array = X.array[,,index]
  Swap = Swap[index]

  aveeff = apply(X.array, 3, function(x) ave.eff(x, blk.proj, Rep))
  index = which(aveeff == max(aveeff, na.rm = TRUE))
  orig.X.array = orig.X.array[,,index]
  X.array= X.array[,,index]
  Swap = Swap[index]

  if(length(dim(X.array)) == 3){
    baleff = apply(X.array, 3, function(x) bal.eff(x, blk.proj, Rep))
    index = which(baleff == min(baleff, na.rm = TRUE))
    orig.X.array = orig.X.array[,,index]
    X.array= X.array[,,index]
    Swap = Swap[index]
  }
  
  return(list(X.trt = X.array, Swap = Swap, orig = orig.X.array))
}




reSwap =
function(des.mat, swap){

  for(i in 1:nrow(swap)){

    temp1 = des.mat[swap[i,1],]
    temp2 = des.mat[swap[i,2],]

    des.mat[swap[i,2],] = temp1
    des.mat[swap[i,1],] = temp2

  }
  
  return(des.mat)
}


reSwap(X.trt, XX$swap[-nrow(XX$swap),])








