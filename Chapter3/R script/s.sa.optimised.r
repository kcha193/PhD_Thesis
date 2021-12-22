    /*
    Rprintf("compareAveEff %f", compareAveEff);
    Rprintf("aveEff %f", aveEff);
    Rprintf("\\n");
    
    Rprintf("compareAveEff %f", min(compare));
    Rprintf("\\n");
     */


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

  arma::umat checkSwap(double aveEff, double oldV, int oldLen,  double rep,
                arma::umat newX, arma::mat Y, arma::mat C, double temp){

    arma::umat newZ = arma::umat(1,1).zeros();

    arma::vec compare = arma::vec(2).ones();
    arma::vec compare1 = arma::vec(2).ones();
    
    arma::mat infoMat = arma::mat(C * arma::trans(newX) * Y * newX * C);

    arma::vec eigval = arma::eig_sym(infoMat);
  
    if(sum(eigval) < 0.0000000001) return(newZ);

    arma::vec conEff = eigval.elem(find(eigval > 0.000001))/rep;
    
    int len = conEff.n_elem;
    
    if(len < oldLen) return(newZ);
    
    double compareAveEff= 1/mean(1/conEff);
    double newV = var(conEff);  
    compare(1) = exp((compareAveEff - aveEff)/temp);
    compare1(1) = exp((oldV - newV)/temp);
           
    if((min(compare) * min(compare1))> 0.25){
      return(newX);
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
  
  arma::mat infoMat;
  arma::umat newZ, newX, bestX;
  arma::vec eigval, conEff, storeAveEff, storeBest;
  arma::urowvec z, oldZ, check;
  double aveEff, oldV, newAveEff, newV, bestAveEff, bestV, theta, delta;
  int oldLen, newLen, bestLen, start;
  
  bestAveEff = 0.0;
  oldZ = arma::urowvec(2).zeros();

  infoMat = arma::mat(C * arma::trans(X) * Y * X * C);
  eigval = arma::eig_sym(infoMat);
  conEff = eigval.elem(find(eigval > 0.000001))/rep;
  oldLen = conEff.n_elem;
  aveEff = 1/(mean(1/conEff));
  oldV = var(conEff);
  
  storeAveEff = arma::vec(tol).zeros();
  storeBest  = arma::vec(tol).zeros();
  bestAveEff = aveEff;
  bestV = oldV;
  bestLen = oldLen;  
  bestX = X;
  delta = 0;
  start = 0;

  for(int i = 0; i < tol; i++){

    newZ = shuffle(Z); 
    z = newZ.row(0);

    if(sum(oldZ == z) == newZ.n_cols) continue;

    oldZ = z;                   

    temp = temp / (1+ cool);

    newX = swap(z[0], z[1], X);  //perform swapping

    if(newX.n_rows == 1) continue;

    newX = checkSwap(aveEff, oldV, oldLen, rep, newX, Y, C, temp);


    if(newX.n_rows != 1){

      infoMat = arma::mat(C * arma::trans(newX) * Y * newX * C);

      eigval = arma::eig_sym(infoMat);
      conEff = eigval.elem(find(eigval > 0.000001))/rep;
      newLen = conEff.n_elem;
      newAveEff = 1/(mean(1/conEff));
      newV = var(conEff); 

      storeAveEff(i) = newAveEff;
      
      if(newAveEff > bestAveEff){
        delta = delta + (newAveEff - bestAveEff);
        start++;

        bestX =  newX;
        bestAveEff = newAveEff; 
        bestV = newV;
        bestLen = newLen;
        storeBest(i) = newAveEff;


        } else if(newAveEff == bestAveEff && newV < bestV){
        theta = (newAveEff - bestAveEff)/(double)(i - start);
        start = i;
        bestX =  newX;
        bestAveEff = newAveEff; 
        bestV = newV;
        bestLen = newLen;
      }
      
        aveEff = newAveEff;
        oldV   = newV;
        oldLen = newLen;
        X = newX;
      }
       //Rcpp::Rcout << "i: " << i << std::endl;
      //Rcpp::Rcout << "temp: " << temp << std::endl;
      //if(temp < 0.0000000000000001) break;

    } 

 
   return Rcpp::List::create(
        Rcpp::Named("X.trt") = bestX,
        Rcpp::Named("nCan")= bestLen,
        Rcpp::Named("maxAveEff")= bestAveEff,
        Rcpp::Named("vCan")= bestV,
        Rcpp::Named("delta")= delta,
        Rcpp::Named("start")= start,
        Rcpp::Named("storeAveEff")=storeAveEff);
  '
  
cpp.test <- cxxfunction(
 signature(X_ = "matrix", Y_ = "matrix", Z_ = "matrix", C_ = "matrix", 
 rep_ = "double", tol_ = "int", temp_ = "double", cool_ = "double"),
  body = code, incl = incltxt, plugin = "RcppArmadillo")

states = lapply(1:nrow(X.trt), function(x) t(cbind(x, (x+1):nrow(X.trt))))
states = states[-length(states)]

#combine all the states into a matrix
states= t(matrix(c(states, recursive=TRUE), nrow = 2))

cpp.states = states - 1
C.trt.mat = diag(ncol(X.trt))

X.trt = X.trt[sample(nrow(X.trt)),]

test(X.trt, blk.proj, Rep = Rep)
#Chain 0
xx =  cpp.test(X.trt, blk.proj, cpp.states, C.trt.mat, Rep, 10000, 10^9, 0)

temp = xx$storeAveEff[-which(xx$storeAveEff==0)]

temp = diff(temp)

temp = temp[which(temp<0)]

(y2 = mean(temp)/log(0.5))
test(X.trt = xx$X.trt, blk.proj, Rep = Rep)

#chain 1
xx =  cpp.test(X.trt, blk.proj, cpp.states, C.trt.mat, Rep, 10000, 10^(-9), 0)

test(X.trt, blk.proj, Rep = Rep)

#chain 2 and onward

xx =  cpp.test(X.trt, blk.proj, cpp.states, C.trt.mat, Rep, 10000, y2, 0)
test(X.trt = xx$X.trt, blk.proj, Rep = Rep)

(y2 = y2/(1+0.25))


colnames(xx$X.trt) = sort(levels(interaction(LETTERS[1:nTrt])))

trt.des = apply(xx$X.trt, 1,  function(x)
colnames(xx$X.trt)[which(as.logical(x))])
 trt.des = as.data.frame(t(sapply(strsplit(trt.des, "\\."), function(x) x)))

design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

design.df = cbind(design.df,
            Ani = t(trt.des),
            Trt = phase1DesignEX1[match(as.character(t(trt.des)),
                          as.character(phase1DesignEX1$Ani)),]$Trt)

design.df

summary.aov.onePhase(design.df,  blk.str = "Run + Tag" , trt.str = "Ani")

summary.aov.twoPhase(design.df,  blk.str2 = "Run" , blk.str1 = "Ani", trt.str = "Trt + Tag")

summary.aov.twoPhase(design.df,  blk.str2 = "Run + Tag" , blk.str1 = "Ani", trt.str = "Trt")

