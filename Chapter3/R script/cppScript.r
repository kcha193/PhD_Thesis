
install.packages(c("Rcpp", "inline", "RcppArmadillo"))

library(Rcpp)
library(inline)
library(RcppArmadillo)
library(compiler)
library(gtools)
library(RcppEigen)

 code <- '
   arma::mat coeff = Rcpp::as<arma::mat>(a);
   arma::mat errors = Rcpp::as<arma::mat>(e);
   int m = errors.n_rows; 
   int n = errors.n_cols;
   arma::mat simdata(m,n);
   simdata.row(0) = arma::zeros<arma::mat>(1,n);
   for (int row=1; row<m; row++) {
     simdata.row(row) = simdata.row(row-1)*trans(coeff)+errors.row(row);
   }
   return Rcpp::wrap(simdata);
   '
  
code<- '
// copy the data to armadillo structures
arma::mat X = Rcpp::as<arma::mat>(X_);
arma::mat Y = Rcpp::as<arma::mat>(Y_);

// calculate the result
arma::mat result = arma::mat(
arma::trans(X) * Y * X
);

double tr =  arma::trace(result);

arma::mat result1 = square(result);
double traceSq =  arma::trace(result1);


arma::vec eigval;
arma::mat eigvec;

arma::eig_sym(eigval, eigvec, result); 

return Rcpp::List::create( 
    Rcpp::Named("info.mat") = result,
    Rcpp::Named("trace")    = tr,
    Rcpp::Named("trace.sq") = traceSq,
    Rcpp::Named("eigval")   = eigval) ;
'


rcppX <- cxxfunction(
 signature(X_ = "matrix", Y_ = "matrix"), code, plugin = "RcppArmadillo" )

code<- '
// copy the data to armadillo structures
arma::mat X = Rcpp::as<arma::mat>(X_);
arma::mat Y = Rcpp::as<arma::mat>(Y_);
int iter = as<int>(iter_);
double rep = as<double>(rep_);

// calculate the result
arma::mat infoMat = arma::mat(arma::trans(X) * Y * X);

//double tr =  arma::trace(infoMat);

//arma::mat squInfoMat = square(infoMat);
//double trSq = arma::trace(squInfoMat);


arma::vec eigval;
arma::mat eigvec;

arma::eig_sym(eigval, eigvec, infoMat); 

double sumEigval = 0.0;
int len = 0;
int length = eigval.n_elem;
for(int i = 0; i < length; i++){
  
  if(eigval[i] < 0.0001) continue;
  
  sumEigval = sumEigval + eigval[i];
  len = len + 1;
}

double aveEff = 1/(rep/(sumEigval/len));

for(int i = 0; i<iter; i++){

  arma::mat newX = shuffle(X);
  
  arma::mat newInfoMat = arma::mat( arma::trans(newX) * Y * newX );

  arma::eig_sym(eigval, eigvec, newInfoMat); 

  double sumEigval = 0.0;
  int len = 0;
  int length = eigval.n_elem;
  for(int i = 0; i < length; i++){
    
    if(eigval[i] < 0.0001) continue;
    
    sumEigval = sumEigval + eigval[i];
    len = len + 1;
  }
  
  double newAveEff = 1/(rep/(sumEigval/len));


  if(newAveEff > aveEff){
    aveEff = newAveEff;
    
    X = newX;     
  }             
}               
 
return Rcpp::List::create( 
    Rcpp::Named("X.trt") = X,
    Rcpp::Named("aveEff")= aveEff);
'

cpp.r.optimised <- cxxfunction(
 signature(X_ = "matrix", Y_ = "matrix", rep_ = "double", iter_ = "int"),
  code, plugin = "RcppArmadillo")





x = diag(10)

cpp.r.optimised(x,x,x, 1, 2)



code<- '
// copy the data to armadillo structures
arma::mat X = Rcpp::as<arma::mat>(X_);

arma::mat Y = shuffle(X);

return Rcpp::wrap(Y);
'

cpp.r.optimised <- cxxfunction(
 signature(X_ = "matrix"),
  code, plugin = "RcppArmadillo")

myfun = function(x) for( i in 1:1000000)  {info.mat = t(x) %*% x %*% x; sum(diag(info.mat))}


myCfun = function(x) for( i in 1:1000000)  {xx = rcppX(x, x, x); sum(xx$tr)}

Rprof("example.out")
myfun(x)
Rprof(NULL)
summaryRprof("example.out") 

Rprof("example.out")
myCfun(x)
Rprof(NULL)
summaryRprof("example.out") 

com.fun = cmpfun(myfun)
Rprof("example.out")
com.fun(x)
Rprof(NULL)
summaryRprof("example.out") 




fx <- cxxfunction(signature( x = "numeric" ),
 'NumericVector xx(x);
return wrap( std::accumulate( xx.begin(), xx.end(), 0.0));',
 plugin = "Rcpp")
 res <- fx( seq( 1, 10, by = 0.5 ) )
 res
 
 
 
 
 
 
code<- '
typedef Eigen::Map<Eigen::MatrixXi> MapMati;

const MapMati B(as<MapMati>(BB));
const MapMati C(as<MapMati>(CC));
return Rcpp::List::create(
  Named("B %*% C") = B * C,
  Named("crossprod(B, C)") = B.adjoint() * C);
'


fprod <-cxxfunction(signature(BB = "MapMati", CC = "MapMati" ), code, plugin="RcppEigen")

B = matrix(rnorm(10), 5,5)
C = matrix(rnorm(10), 5,5)

prodCpp(B,C)

A <- matrix(1:6, ncol=2)

 ftrans(A)
 
Rsamp <- function(X) {
  stopifnot(is.numeric(X <- as.matrix(X)),
            (nc <- ncol(X)) > 1L,
            all(X >= 0))
  apply(X, 1L, function(x) sample(nc, size=1L, replace=FALSE, prob=x+1e-10))
}


code<-'
  typedef Eigen::ArrayXd   Ar1;
  typedef Eigen::ArrayXXd  Ar2;
  typedef Eigen::Map<Ar2> MAr2;

  const MAr2 X(as<MAr2>(X_));
  int m(X.rows()), n(X.cols()), nm1(n - 1);
  Ar1     samp(m);
  RNGScope sc; // Handle GetRNGstate()/PutRNGstate()
  for (int i=0; i < m; ++i) {
    Ar1 ri(X.row(i));
    std::partial_sum(ri.data(), ri.data() + n, ri.data());
    ri /= ri[nm1];    // normalize to sum to 1
    samp[i] = (ri < ::unif_rand()).count() + 1;
  }
  return Rcpp::wrap(samp);
'

ftrans <-cxxfunction(signature(X_="MAr2"), code, plugin="RcppEigen")


RcppSamp <- function(X) {
  stopifnot(is.numeric(X <- as.matrix(X)),
            (nc <- ncol(X)) > 1L,
            all(X >= 0))
  ftrans(X)
}

set.seed(1234321)
X <- matrix(runif(100000 * 5), ncol=5)
benchmark(Rsamp(X), RcppSamp(X), replications=5,
          columns=c("test", "elapsed", "relative", "user.self"))

