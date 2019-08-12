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
// [[Rcpp::export]]
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

//repeat the same value
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

//compute and sum over quadratic form (y_{i} - mu)B(y_{i} - mu)'
// [[Rcpp::export]]
double sumQFCentSq(arma::mat y, arma::vec mu, arma::mat B) {
  int n = y.n_cols;
  arma::vec sum = arma::zeros(1);
  for (unsigned i = 0; i < n; i++) {
    sum = sum + (y.col(i) - mu).t()*B*(y.col(i) - mu);
  }
  return(sum(0));
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

//convert x values (locations of points on contours) and center point to
//y values (distances from center points)
//[[Rcpp::export]]
arma::mat xToY(arma::cube x, arma::vec C) {
  int nPts = x.n_cols;
  int nSamp = x.n_slices;
  arma::mat y = arma::zeros(nPts, nSamp);
  for (unsigned i = 0; i < nSamp; i++) {
    arma::mat xCurr = x.slice(i);
    for (unsigned j = 0; j < nPts; j++) {
      y(j, i) = ptDist(xCurr.col(j), C);
    }
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
      if (ea(1) == pt(1)) { //test point and edge must be at same height for point to be in
        if (pt(0) >= std::min(ea(0), eb(0)) & pt(0) <= std::max(ea(0), eb(0))) {
          return(TRUE);
        }
      }

    //vertical edge
    } else if (ea(0) == eb(0)) {
      if (ea(0) > pt(0)) { //line must be to the right
        if ((pt(1) >= std::min(ea(1), eb(1))) &&
            (pt(1) <= std::max(ea(1), eb(1)))) {//test point y must be between edge end points y's
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
      //find x-value at which the edge intersects the height of the test
      //ray. Then check if that x-value is on or to the right of the point
      double m = (eb(1) - ea(1))/(eb(0) - ea(0));
      double b = eb(1) - m*eb(0);
      double xInter = (pt(1) - b)/m;
      if (abs(xInter - pt(0)) <= 1e-10) { //test point is (approximately) on edge
        return(TRUE);
      } else if ((xInter > pt(0)) && (xInter >= std::min(ea(0), eb(0)))
            && (xInter <= std::max(ea(0), eb(0)))) {
        //x intersection needs to be on edge line segment and to the right of test pt
        count += 1;
      }
    }
  }
  
  if (count % 2 == 0) {
    return(FALSE);
  } else {
    return(TRUE);
  }
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
List RunMCMC(int n_iter, int w, arma::cube x, arma::vec C,
             arma::vec mu0,  arma::mat lambda0, arma::mat S0, int nu0,
             arma::mat Sigma,
             arma::vec mu0C,  arma::mat lambda0C, arma::mat S0C, int nu0C,
             arma::mat SigmaC, 
             double CxSD, double CySD, arma::mat kernHat) {

  Rcout << "starting " << 0 << std::endl;
  
  //constants
  int nVecs = x.n_cols;
  int nObs = x.n_slices;
  int nuN = nu0 + nObs;

  //storage vectors and matrices
  arma::mat muStore(nVecs, n_iter/w);
  arma::cube sigmaStore(nVecs, nVecs, n_iter/w);

  //Cx and Cy acceptance rates
  double accRateCx = 0;
  double accRateCy = 0;

  //inverses
  arma::mat SigmaInv = inv(Sigma);
  arma::mat lambda0Inv = inv(lambda0);

  //compute initial y and muC
  arma::mat y = xToY(x, C);
  arma::vec muC = C;
  
  //create prop vectors
  arma::vec CProp = C;

  
  //Rcout << "prelim done " << 0 << std::endl;
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
    //Rcout << "gibbs for mu done" << j << std::endl;
  
   //////////Gibbs for Sigma/////////////
   arma::mat Sn = S0 + sumQFCent(y, mu);
   arma::mat Sigma = riwish(nuN, inv(Sn));
   if (j%w == 0) {
     sigmaStore.slice(j/w) = Sigma;
   }
  //Rcout << "gibbs for sigma done" << j << std::endl;
   
  ////////////Metropolis for Cx///////
  arma::vec Cx = arma::randn(1)*CxSD + C(0);
  CProp(0) = Cx(0);
  if (ptInPoly(kernHat, CProp)) {
    //Rcout << "Cx after ptInpoly" << j << std::endl;
    arma::mat yxProp = xToY(x, CProp);
    bool acceptCx = false;
    double logRCx = (- 0.5*sumQFCentSq(yxProp, mu, inv(Sigma))
                     - 0.5*sumQFSq(CProp - muC, inv(SigmaC))
                     + 0.5*sumQFCentSq(y, mu, inv(Sigma))
                     + 0.5*sumQFSq(C - muC, inv(SigmaC)));
    //Rcout << "have logRCx" << logRCx << std::endl;
    acceptCx = logAccept(logRCx);
    //Rcout << "acceptCx" << acceptCx << std::endl;
    
    if (acceptCx) {
      //Rcout << "in accept Cx" << j << std::endl;
      C = CProp;
      accRateCx++;
      y = yxProp;
    } else {
     //Rcout << "in reject Cx" << j << std::endl;
      CProp = C;
    }
  } else {
    //Rcout << "rejected bc ptInPoly Cx" << j << std::endl;
    CProp = C;
  }

  ////////////Metropolis for Cy///////
  arma::vec Cy = arma::randn(1)*CySD + C(1);

  CProp(1) = Cy(0);
    if (ptInPoly(kernHat, CProp)) {
    arma::mat yyProp = xToY(x, CProp);
    bool acceptCy = false;
    double logRCy = (- 0.5*sumQFCentSq(yyProp, mu, inv(Sigma))
                     - 0.5*sumQFSq(CProp - muC, inv(SigmaC))
                     + 0.5*sumQFCentSq(y, mu, inv(Sigma))
                     + 0.5*sumQFSq(C - muC, inv(SigmaC)));
    acceptCy = logAccept(logRCy);
    //Rcout << "acceptCy" << acceptCy << std::endl;
    
    if (acceptCy) {
      //Rcout << "in accept Cy" << j << std::endl;
      C = CProp;
      accRateCy++;
      //Rcout << "accRateCy" << accRateCy << std::endl;
      
      y = yyProp;
    } else {
      //Rcout << "in reject Cy" << j << std::endl;
      CProp = C;
    }
  } else {
    //Rcout << "rejected bc ptInPoly" << j << std::endl;
    CProp = C;
  }
  //Rcout << "metropolis for Cy done" << j << std::endl;
  

  ////////Gibbs for muC/////////////


  //////////Gibbs for SigmaC/////////////



  }

  //finalize acceptance rate
  accRateCx = accRateCx/n_iter;
  accRateCy = accRateCy/n_iter;
  
  
  //return values
  List res;
  res["mu"] = muStore; res["Sigma"] = sigmaStore; 
  res["accRateCx"] = accRateCx; res["accRateCy"] = accRateCy;
  return(res);
}

