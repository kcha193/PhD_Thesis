
code<- '
// copy the data to armadillo structures
arma::mat X = Rcpp::as<arma::mat>(X_);
arma::mat Y = Rcpp::as<arma::mat>(Y_);
int iter = as<int>(iter_);
double rep = as<double>(rep_);

//shuffle the rows of the treatment design matrix

// calculate the information matrix
arma::mat infoMat = arma::mat(arma::trans(X) * Y * X);

arma::mat eigvec;

arma::vec eigval = arma::eig_sym(infoMat);
arma::rowvec compare1 = X.row(0);
arma::rowvec compare2 = X.row(4);

eigval = eigval.elem(find(eigval > 0.000001) );
return Rcpp::wrap(sum(compare1==compare2));
'
cpp.test <- cxxfunction(
 signature(X_ = "matrix", Y_ = "matrix", rep_ = "double", iter_ = "int"),
  code, plugin = "RcppArmadillo")

cpp.test(X.trt, blk.proj, Rep, 10000000)

#ranomdised method in c++ implementation
code<- '
// copy the data to armadillo structures
arma::umat X = Rcpp::as<arma::umat>(X_);
arma::mat Y = Rcpp::as<arma::mat>(Y_);
int iter = as<int>(iter_);
double rep = as<double>(rep_);

//Declare variables
arma::mat infoMat, newInfoMat, squInfoMat, squNewInfoMat;
arma::umat newX;
arma::vec conEff, eigval, newConEff;
double aveEff, newAveEff, newTr, newTrSq, newV, oldV,  tr, trSq;
double sumEigval = 0.0;
int newLen;
int oldLen = 0;
       
//shuffle the rows of the treatment design matrix
X = shuffle(X);

// calculate the information matrix
infoMat = arma::mat(arma::trans(X) * Y * X);

int counter = 0;

while(counter < iter){
  
  newX = shuffle(X);

  counter++;
  
  newInfoMat = arma::mat( arma::trans(newX) * Y * newX );

  //calculate the trace of the matrix
 tr =  arma::trace(infoMat);

 newTr =  arma::trace(newInfoMat);

  //compare the trace of the matrix
  if(newTr>tr){
  	tr = newTr;
  	X = newX;

  } else if (newTr == tr){

	eigval = arma::eig_sym(infoMat);

	conEff = eigval.elem(find(eigval > 0.000001) );
	
	sumEigval = sum(conEff);

  oldLen = conEff.n_elem;
  
	 aveEff = 1/(rep/(sumEigval/oldLen));


	 eigval = arma::eig_sym(newInfoMat); 

	  newConEff = eigval.elem(find(eigval > 0.000001) );

    sumEigval = sum(newConEff);

    newLen = newConEff.n_elem;
    
	  //check for the number of non-zero eigenvalues
	  if(newLen < oldLen) continue;

	  //calculate the average efficiency factor
	   newAveEff = 1/(rep/(sumEigval/newLen));

	  //compare the average efficiency factors
	  if(newAveEff > aveEff){
  		oldLen = newLen;
  		aveEff = newAveEff;
  		conEff = newConEff;
  
  		X = newX;

	  } else if( newAveEff == aveEff){
		//calculate the trace of the sqaured matrix
		squInfoMat = square(infoMat);
		trSq = arma::trace(squInfoMat);

		squNewInfoMat = square(newInfoMat);
		newTrSq = arma::trace(squNewInfoMat);

		if(newTrSq < trSq){
			trSq = newTrSq;
			oldLen = newLen;
			aveEff = newAveEff;
			conEff = newConEff;

			X = newX;

		} else if(newTrSq == trSq){


		 oldV = var(conEff);
		 newV = var(newConEff);

		if(oldV < newV){
		  oldLen = newLen;
		  aveEff = newAveEff;
		  conEff = newConEff;

		  X = newX;

		}
	  }

	}
}
}

return Rcpp::List::create(
    Rcpp::Named("X.trt") = X,
    Rcpp::Named("aveEff")= aveEff,
    Rcpp::Named("conEff")= conEff/rep);
'

cpp.r.optimised <- cxxfunction(
 signature(X_ = "matrix", Y_ = "matrix", rep_ = "double", iter_ = "int"),
  code, plugin = "RcppArmadillo")


r.optimised =
function(X.trt, blk.proj, nItr = 10000000){
  cat("Finding the optimal design using random sampling:\n")

  X.trt= X.trt[sample(nrow(X.trt)),]
  info.mat = t(X.trt) %*% blk.proj %*%X.trt
   
  temp.trace = tr(info.mat)
  temp.squ.trace = tr(info.mat^2)

  e.va = eigen(info.mat)$va
  ave.e.va = (length(e.va[-which(e.va<1e-7)])/sum(1/e.va[-which(e.va<1e-7)]))/Rep
  sd.e.va = sd(e.va[-which(e.va<1e-7)])

  if(is.na(sd.e.va))  sd.e.va = 0

  pb <- txtProgressBar(min = 0, max = nItr, style = 3)
  for( i in 1:nItr){
    setTxtProgressBar(pb, i)

    X.trt1 = X.trt[sample(nrow(X.trt)),]
    info.mat = t(X.trt1) %*% blk.proj %*%X.trt1
 
    new.temp.trace = tr(info.mat)

    if( new.temp.trace>temp.trace) {
      temp.trace = new.temp.trace
      X.trt = X.trt1
    } else if(new.temp.trace == temp.trace){

      e.va = eigen(info.mat)$va
      temp.ave.e.va = (length(e.va[-which(e.va<1e-7)])/sum(1/e.va[-which(e.va<1e-7)]))/Rep

      if(temp.ave.e.va > ave.e.va){
        ave.e.va = temp.ave.e.va

        X.trt = X.trt1

      } else if(temp.ave.e.va == ave.e.va){

        new.temp.squ.trace = tr(info.mat^2)

        if(new.temp.squ.trace < temp.squ.trace){
          temp.squ.trace = new.temp.squ.trace
          X.trt = X.trt1
        }else if(new.temp.squ.trace == temp.squ.trace) {
          temp.sd.e.va = sd(e.va[-which(e.va<1e-7)])

            if(is.na(temp.sd.e.va))  temp.sd.e.va = 0

          if(temp.sd.e.va < sd.e.va){
            sd.e.va = temp.sd.e.va
            X.trt = X.trt1
          }

        }
      }
    }

  }
  close(pb)

  return(X.trt1)
}


cpp.r.optimised =
function(X.trt, blk.proj, nItr = 10000000){
  #cat("Finding the optimal design using random sampling:\n")

  X.trt= X.trt[sample(nrow(X.trt)),]
  temp = rcppX( X.trt, blk.proj)
  info.mat = temp$info.mat
   
  temp.trace = temp$trace   
  temp.squ.trace = temp$trace.sq 
  
  e.va = as.numeric(temp$eigval)
  e.va = e.va[-which(e.va<1e-7)]
  ave.e.va = 1/mean(Rep/e.va)
  sd.e.va = sd(e.va)

  if(is.na(sd.e.va))  sd.e.va = 0

  #pb <- txtProgressBar(min = 0, max = nItr, style = 3)
  for( i in 1:nItr){
   # setTxtProgressBar(pb, i)

    X.trt1 = X.trt[sample(nrow(X.trt)),]
    temp = rcppX(X.trt1, blk.proj)
    info.mat = temp$info.mat
   
   new.temp.trace = temp$trace   

    if( new.temp.trace>temp.trace) {
      temp.trace = new.temp.trace
      X.trt = X.trt1
    } else if(new.temp.trace == temp.trace){

      e.va = as.numeric(temp$eigval)
      e.va = e.va[-which(e.va<1e-7)]
      temp.ave.e.va = 1/mean(Rep/e.va)

      if(temp.ave.e.va > ave.e.va){
        ave.e.va = temp.ave.e.va

        X.trt = X.trt1

      } else if(temp.ave.e.va == ave.e.va){

        new.temp.squ.trace = temp$trace.sq 
   
        if(new.temp.squ.trace < temp.squ.trace){
          temp.squ.trace = new.temp.squ.trace
          X.trt = X.trt1
        }else if(new.temp.squ.trace == temp.squ.trace) {
          temp.sd.e.va = sd(e.va)

            if(is.na(temp.sd.e.va))  temp.sd.e.va = 0

          if(temp.sd.e.va < sd.e.va){
            sd.e.va = temp.sd.e.va
            X.trt = X.trt1
          }

        }
      }
    }

  }
 # close(pb)

  return(X.trt1)
}

