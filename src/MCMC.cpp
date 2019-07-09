#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// draw from multivariate normal dist using Cholesky decomp. method
arma::mat mvrnormCpp(arma::vec mu, arma::mat Sigma) {
  int nCol = Sigma.n_cols;
  arma::mat u = arma::randn(1, nCol);
  arma::mat temp = mu.t() + u*arma::chol(Sigma);
  return temp.t();
}

// compute whether to accept or reject proposal based on
// log of acceptance rate
bool logAccept(double logR) {
  if (logR > 0.0) { //accept with prob min(r, 1)
    logR = 0.0;
  }
  double logU = log(R::runif(0, 1));
  if (logU < logR) {
    return(true);
  }
  else {
    return(false);
  }
}

// template to compute mean
template <class T>
inline double GetMean(T& x ) {
  return mean(x);
}

// compute row means
arma::vec RowMean(arma::mat m) {
  int nRow = m.n_rows;
  arma::vec out(nRow);
  for (int i = 0; i < nRow; i++) {
    arma::rowvec c = m.row(i);
    out(i) = GetMean(c);
  }
  return(out);
}

// compute and sum over quadratic form a_{i}'ba_{i}
double sumQFSq(arma::mat a, arma::mat b) {
  arma::vec sum = arma::zeros(1);
  int n = a.n_cols;
  for (unsigned i = 0; i < n; i++) {
    sum = sum + a.col(i).t()*b*a.col(i);
  }
  return(sum(0));
}

// compute and sum over quadratic form a_{i}'bc_{i}
double sumQFCon(arma::mat a, arma::mat b, arma::colvec c) {
  arma::vec sum = arma::zeros(1);
  int n = a.n_cols;
  for (unsigned i = 0; i < n; i++) {
    sum = sum + a.col(i).t()*b*c;
  }
  return(sum(0));
}

//repeate the same value
arma::colvec rep(double a, int n) {
  arma::colvec A = arma::zeros(n);
  for (int i = 0; i < n; i++) {
    A(i) = a;
  }
  return(A);
}

//compute and sum over quadratic form (y_{i} - mu)(y_{i} - mu)'
// [[Rcpp::export]]
arma::mat sumQFCent(arma::mat y, arma::vec mu) {
  int n = y.n_rows;
  int m = y.n_cols;
  arma::mat sum = arma::zeros<arma::mat>(n,n);
  for (unsigned i = 0; i < m; i++) {
    sum = sum + (y.col(i) - mu)*(y.col(i) - mu).t();
  }
  return(sum);
}

//inverse wishart samples, computed as described in Hoff pg 110
// [[Rcpp::export]]
arma::mat riwish(int nu, arma::mat SInv) {
  arma::mat SigmaInv;
  int n = SInv.n_cols;
  arma::vec mu_zero = arma::zeros(n);
  for (unsigned i = 0; i < nu; i++) {
    arma::vec zi = mvrnormCpp(mu_zero, SInv);
    if (i == 0) {
      SigmaInv  = zi*zi.t();
    } else {
      SigmaInv = SigmaInv + zi*zi.t();
    }
  }
  arma::mat Sigma = inv(SigmaInv);

  return(Sigma);
}


//' Run MCMC to Fit Contour Model
//' @param n_iter   number of iterations to run the MCMC
//' @param u a matrix of observed lengths of dimension number
//'          of vectors by number of observed contours
//' @param mu vector of the same length as the number of lines which specifies
//'           the values from which each element of \code{mu} will be initialized
//'           in the MCMC.
//' @param mu0 vector of the same length as the number of lines which specifies
//'            the prior mean for \code{mu}.
//' @param lambda0 matrix of the same dimension as the number of lines which
//'                specifices the prior covariance matrix for \code{mu}.
//' @param S0 matrix BLAH
//'
//' @param w Integer specifying how many samples of the parameters will be
//'           maintained. Samples from every wth iteration is stored.
//'
//' @return List of length X that gives BLAH
// [[Rcpp::export]]
List RunMCMC(int n_iter,
             arma::mat y, arma::vec mu0,  arma::mat lambda0,
             arma::mat S0, int nu0,
             arma::mat Sigma_ini, int w) {
  //constants
  int nVecs = y.n_rows;
  int nObs = y.n_cols;
  int nu = nu0 + nVecs;

  //storage vectors and matrices
  arma::mat muStore(nVecs, n_iter/w);
  arma::cube sigmaStore(nVecs, nVecs, n_iter/w);

  //inverses
  arma::mat Sigma = Sigma_ini;
  arma::mat SigmaInv = inv(Sigma);
  arma::mat lambda0Inv = inv(lambda0);
  
  arma::vec mu = mu0; //testing

  // // /////////main loop/////////////
  for (unsigned j = 0; j < n_iter; j++) {

    ////////Gibbs for mu/////////////
    arma::vec yMean = RowMean(y);
    arma::mat lambdaN = inv(lambda0Inv + nObs*SigmaInv);
    arma::vec muN = (inv(lambda0Inv + nObs*SigmaInv)*
      (lambda0Inv*mu0 + nObs*SigmaInv*yMean));
    arma::mat mu = mvrnormCpp(muN, lambdaN);
    if (j%w == 0) {
      muStore.col(j/w) = mu;
    }

   //////////Gibbs for Sigma/////////////
   arma::mat Sn = S0 + sumQFCent(y, mu);
   arma::mat Sigma = riwish(nu, inv(Sn));
   if (j%w == 0) {
     sigmaStore.slice(j/w) = Sigma;
   }
  }

  //return values
  List res;
  res["mu"] = muStore; res["Sigma"] = sigmaStore;
  return(res);
}

