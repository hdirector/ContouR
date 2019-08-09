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

// find the distance between two points
// [[Rcpp::export]]
double ptDist(arma::vec p1, arma::vec p2) {
  return(sqrt(arma::sum(arma::square(p1 - p2))));
}

// convert x values (locations of points on contours) and center point to
// y values (distances from center points)
// [[Rcpp::export]]
arma::vec xToY(arma::mat x, arma::vec C) {
  int n = x.n_cols;
  arma::vec y =  arma::zeros(n);
  for (unsigned i = 0; i < n; i++) {
    y(i) = ptDist(x.col(i), C);
  }
  return(y);
}


//check if a point is in a polygon using the ray-casting approach
//e.g., https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
// [[Rcpp::export]]
bool ptInPoly(arma::mat poly, arma::vec pt) {
  int n = poly.n_cols;
  int count = 0; //count of intersection with test ray and polygon edges

  //main section: loop over the edges and count intersections with test ray

  
  for (unsigned i = 0; i < n - 1; i++) {
    arma::vec ea = poly.col(i);
    arma::vec eb = poly.col(i + 1);
    //horizontal edge
    if (ea(1) == eb(1)) {  
      Rcout << "horizontal edge case:" << i << std::endl;
      if (ea(1) == pt(1)) { //test point and edge must be at same height for point to be in
        if (pt(0) >= std::min(ea(0), eb(0)) & pt(0) <= std::max(ea(0), eb(0))) {
          return(TRUE);
        }
      }
    
    //vertical edge  
    } else if (ea(0) == eb(0)) { 
      Rcout << "vertical edge case:" << i << std::endl;
      if (ea(0) > pt(0)) { //line must be to the right
        Rcout << "line to right:" << 22 << std::endl;
        if ((pt(1) >= std::min(ea(1), eb(1))) &&
            (pt(1) <= std::max(ea(1), eb(1)))) {//test point y must be between edge end points y's
          Rcout << "in count:" << 22 << std::endl;
          count += 1;
        }
      } else if (ea(0) == pt(0)) {
        if ((pt(1) >= std::min(ea(1), eb(1))) &&
            (pt(1) <= std::max(ea(1), eb(1)))) { //test point on boundary line
          return(TRUE);
        }
      }
      
    //typical edge  
    } else {  
      Rcout << "typical edge case:" << i << std::endl;
      //find x-value at which the edge intersects the height of the test
      //ray. Then check if that x-value is on or to the right of the point
      double m = (eb(1) - ea(1))/(eb(0) - ea(0));
      double b = eb(1) - m*eb(0);
      double xInter = (pt(1) - b);
      if (xInter >= pt(0)) {
        count += 1;
      }
    }
  }
 
  
  if (count % 2 == 0) {
    Rcout << "count even :" << count << std::endl;
    return(FALSE);
  } else {
    Rcout << "count odd:" << count << std::endl;
    return(TRUE);
  }
}



// [[Rcpp::export]]
double test(int x) {
  if ( x % 2 == 0) {
    Rcout << "The value is even" <<  x<< std::endl;
  } else {
    Rcout << "The value is odd" << x << std::endl;
  }
  return(0);
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
List RunMCMC(int n_iter, int w, arma::mat x, arma::vec C,
             arma::vec mu0,  arma::mat lambda0, arma::mat S0, int nu0,
             arma::mat Sigma_ini,
             arma::vec mu0C,  arma::mat lambda0C, arma::mat S0C, int nu0C,
             arma::mat SigmaIniC) {

  //constants
  int nVecs = x.n_rows;
  int nObs = x.n_cols;
  int nuN = nu0 + nObs;

  //storage vectors and matrices
  arma::mat muStore(nVecs, n_iter/w);
  arma::cube sigmaStore(nVecs, nVecs, n_iter/w);

  //inverses
  arma::mat Sigma = Sigma_ini;
  arma::mat SigmaInv = inv(Sigma);
  arma::mat lambda0Inv = inv(lambda0);

  //compute initial y
  arma::vec y = xToY(x, C);

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
   arma::mat Sigma = riwish(nuN, inv(Sn));
   if (j%w == 0) {
     sigmaStore.slice(j/w) = Sigma;
   }
  }

  //return values
  List res;
  res["mu"] = muStore; res["Sigma"] = sigmaStore;
  return(res);
}

