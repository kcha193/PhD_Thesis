
 code<- '
  // copy the data to armadillo structures
  arma::mat X = Rcpp::as<arma::mat>(X_);
   arma::cube xMats(X.n_rows, X.n_cols, 1);
  xMats.zeros();
  
  xMats.resize(5 + 2, 2, 1);


   
  return Rcpp::wrap(xMats);
  '
cpp.test <- cxxfunction(
 signature(X_ = "matrix"), body = code,  plugin = "RcppArmadillo")
 
cpp.test(X.trt)

incltxt <- '
    arma::umat swap(arma::umat Z, arma::umat X) {

      arma::umat newX = X;
      arma::umat newZ = arma::umat(1,1).zeros();
      int length = Z.n_rows;
      for(int i = 1; i<length; i++){	 
        arma::urowvec z = Z.row(i);
        
        int s1 = z[0];
        int s2 = z[1];    
        
        arma::urowvec X1 = X.row(s1);
        arma::urowvec X2 = X.row(s2);
                
        if(sum(X1 == X2) == X.n_cols) return(newZ);
        
        newX.row(s1) = X.row(s2);
        newX.row(s2) = X.row(s1);
      }
      return(newX);
    }
    
  arma::urowvec checkSwap(double aveEff, double oldV, int oldLen,  double rep, arma::umat curZ,
                arma::urowvec z, arma::umat X, arma::mat Y, arma::mat C, double temp){
    
    int current = curZ.n_rows;

    curZ.insert_rows(current, z);
    arma::umat newX = swap(curZ, X);
    
    arma::urowvec newZ = arma::urowvec(1).zeros();    
    
    if(newX.n_rows == 1) return(newZ);

    arma::vec compare = arma::vec(2).ones();
    arma::vec compare1 = arma::vec(2).ones();
    
    arma::mat infoMat = arma::mat(C * arma::trans(newX) * Y * newX * C);

    arma::vec eigval = arma::eig_sym(infoMat);
  
    arma::vec conEff = eigval.elem(find(eigval > 0.000001))/rep;
    
    int len = conEff.n_elem;
    
    if(len < oldLen) return(newZ);
    
    double compareAveEff= 1/mean(1/conEff);
    double newV = var(conEff);  
    compare(1) = exp((compareAveEff - aveEff)/temp);
    compare1(1) = exp((oldV - newV)/temp);
           
    if((min(compare)*min(compare1)) > 0.95){ 
      return(z);
    } else{
      return(newZ);
    } 
    
  } 
'


code<- '
  // copy the data to armadillo structures
  arma::umat X = Rcpp::as<arma::umat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::umat Z = Rcpp::as<arma::umat>(Z_);
  arma::mat C = Rcpp::as<arma::mat>(C_);
  double rep = as<double>(rep_);
  int tol = as<int>(tol_);
  double temp = as<double>(temp_);
  double cool = as<double>(cool_);
  arma::ucube xMats(1, 2, 1);
  
  arma::mat infoMat;
  arma::umat newZ, curZ, oldZ, newX, bestCurZ, newCurZ;
  arma::vec eigval;
  arma::vec conEff;
  arma::urowvec z, check;
  double aveEff, oldV, newAveEff, newV, bestAveEff, bestV;
  int oldLen, newLen, bestLen;
  
  bestAveEff = 0.0;
  curZ = arma::umat(1, 2);
  curZ.zeros();
  xMats.slice(0) = curZ;
  
  infoMat = arma::mat(C * arma::trans(X) * Y * X * C);
  eigval = arma::eig_sym(infoMat);
  conEff = eigval.elem(find(eigval > 0.000001))/rep;
  oldLen = conEff.n_elem;
  aveEff = 1/(mean(1/conEff));
  oldV = var(conEff);
  
  bestCurZ =  newCurZ;
  bestAveEff = aveEff; 
  bestV = oldV;
  bestLen = oldLen;  
  
  int lenList = 1;
  int counter = 0;
  int length = Z.n_rows; 
  
  arma::ucube xMats1(2, 2, 1);  
  xMats1.zeros();
  xMats.zeros(); 
  int k = 0;
  
  for(k = 0; k < 2; k++){
  
    for( int j = 0; j < lenList; j++){
      temp = temp * cool;
      curZ =  xMats.slice(j);
      for(int i = 0; i < length; i++){ 
        
        z = Z.row(i);            
        
        if(sum(curZ.row(curZ.n_rows-1) == z) == Z.n_cols) continue;
   
        check = checkSwap(aveEff, oldV, oldLen, rep, curZ, z, X, Y, C, temp);
        
        /*
        Rprintf("bestAveEff %f", bestAveEff);
        Rprintf("\\n");
        Rprintf("aveEff %f", aveEff);
        Rprintf("\\n");
        Rprintf("iter %d", j);
        Rprintf("\\n");
  
         Rprintf("iter %d", i);
        Rprintf("\\n");
         Rprintf("iter %d", k);
        Rprintf("\\n");
          Rprintf("check %d", z[0]);
          Rprintf("check %d", z[1]);
  
         Rprintf("\\n");
          Rprintf("check %d", check[0]);
        
          Rprintf("check %d", check[1]);
         Rprintf("\\n");
          */
          
        if(check.n_elem !=1){   
               
          int current = curZ.n_rows;
          oldZ =  curZ;      
          curZ.insert_rows(current, z); 
          xMats1.insert_slices(0, 1);
       
          xMats1.slice(0) =  curZ; 
          counter++; 
        
          newX = swap(curZ, X);
          
          infoMat = arma::mat(C * arma::trans(newX) * Y * newX * C);
          eigval = arma::eig_sym(infoMat);
          conEff = eigval.elem(find(eigval > 0.000001))/rep;
          newLen = conEff.n_elem;
          newAveEff = 1/(mean(1/conEff));
          newV = var(conEff); 
          
          if(newAveEff > bestAveEff){
            bestCurZ =  curZ;
            bestAveEff = newAveEff; 
            bestV = newV;
            bestLen = newLen;   
          } else if(newAveEff == bestAveEff && newV < bestV){
            bestCurZ =  curZ;
            bestAveEff = newAveEff; 
            bestV = newV;
            bestLen = newLen;     
          }
            curZ =  oldZ;
            aveEff = newAveEff;
            oldV   = newV;
            oldLen = newLen;      
        } 
       
       }  
     }
   
    xMats = xMats1;
    xMats1.reset();
    xMats1.resize(k + 3, 2, 1);
    xMats1.zeros();
    lenList = counter;
    counter = 0;
   }
   return Rcpp::List::create(
        Rcpp::Named("X.trt") = counter,
        Rcpp::Named("nCan")= bestLen,
        Rcpp::Named("maxAveEff")= bestAveEff,
        Rcpp::Named("vCan")= bestV,
        Rcpp::Named("swap")= xMats);
  '
  
cpp.test <- cxxfunction(
 signature(X_ = "matrix", Y_ = "matrix", Z_ = "matrix", C_ = "matrix", 
 rep_ = "double", tol_ = "int", temp_ = "double", cool_ = "double"),
  body = code, incl = incltxt, plugin = "RcppArmadillo")

   xx =  cpp.test(X.trt[sample(nrow(X.trt)),], blk.proj, cpp.states, 
   C.trt.mat, Rep, 100, 1000, 0.99)   


 
  states = lapply(1:nrow(X.trt), function(x) t(cbind(x, (x+1):nrow(X.trt))))
  states = states[-length(states)]
  
  #combine all the states into a matrix
  states= t(matrix(c(states, recursive=TRUE), nrow = 2))
       
  cpp.states = states - 1
   C.trt.mat = diag(ncol(X.trt))







