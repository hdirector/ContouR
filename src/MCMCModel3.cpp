#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// draw from multivariate normal dist using Cholesky decomp. method
// [[Rcpp::export]]
arma::vec mvrnormCpp(arma::vec mu, arma::mat Sigma) {
  int nCol = Sigma.n_cols;
  arma::mat u = arma::randn(1, nCol);
  arma::mat temp = mu.t() + u*arma::chol(Sigma);
  arma::vec fin = temp.t();
  return(fin);
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


// find the distance between two points
// [[Rcpp::export]]
double ptDist(arma::vec p1, arma::vec p2) {
  return(sqrt(arma::sum(arma::square(p1 - p2))));
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

//projected point
//[[Rcpp::export]]
arma::vec projDist(arma::vec Cx, arma::vec Cy, arma::vec theta, arma::vec pt) {
  double pi = 3.14159;
  double eps = .00001;
  arma::vec w(1);
  arma::vec y(1);
  //vertical line
  if  (((theta(0) >=  -eps) & (theta(0) <=  eps)) |
       ((theta(0) >= pi - eps) & (theta(0) <= pi + eps))) {
    w = abs(pt(1) - Cy);
    y = abs(pt(0) - Cx);
  //horizontal line
  } else if (((theta(0) >=  pi/2 - eps) & (theta(0) <=  pi/2 + eps)) |
             ((theta(0) >= 3*pi/2 - eps) & (theta(0) <= 3*pi/2 + eps))) {
    arma::vec ptAbs = arma::abs(pt);
    w = abs(pt(0) - Cx);
    y = abs(pt(1) - Cy);
  //typical case
  } else {
    //equation of generating line (line that extends from (Cx, Cy) along theta)
    arma::vec m = sin(theta)/cos(theta);
    arma::vec b = Cy - m*Cx;
    //equation of perpendicular line
    arma::vec mProj = -(1/m);
    arma::vec bProj = pt(1) - mProj*pt(0); 
    //intersection of generating line and perpendicular line
    arma::vec xInter = (b - bProj)/(mProj - m);
    arma::vec yInter = m*xInter + b;
    //w value (perpendicular length)
    w = sqrt(arma::square(pt(0) - xInter) + arma::square(pt(1) - yInter));
    y = sqrt(arma::square(Cx(0) - xInter) + arma::square(Cy(0) - yInter));
  }
  
  arma::vec ret(2);
  ret(0) = w(0);
  ret(1) = y(0);
  return(ret);
}

//Compute W and Y values
//[[Rcpp::export]]
List XToWY(arma::cube x, arma::vec Cx, arma::vec Cy, arma::vec thetas) {
  int p = x.n_cols;
  int nSamp = x.n_slices;
  arma::mat w(p, nSamp);
  arma::mat y(p, nSamp);
  arma::vec temp(2);
  arma::vec ptCurr(2); //TO DO, write this in a less dumb way
  arma::vec thetaCurr(1);
  for (unsigned i = 0; i < p; i++) {
    thetaCurr(0) = thetas(i);
    for (unsigned j = 0; j < nSamp; j ++) {
      ptCurr(0) = x(0,i, j);//TO DO, write this in a less dumb way
      ptCurr(1) = x(1,i, j);//TO DO, write this in a less dumb way
      
      temp = projDist(Cx, Cy, thetaCurr, ptCurr);
      w(i, j) = temp(0);
      y(i, j) = temp(1);
      //Rcout << "inProjDistLoop " << i << std::endl;
      //Rcout << "inProjDistLoop " << j << std::endl;
      
    }
  }
    
  List wy; 
  wy["w"] = w;
  wy["y"] = y;
  return(wy);
}


//Compute Sigma based on sigma and kappa
//[[Rcpp::export]]
arma::mat compSigma(arma::vec sigma, arma::vec kappa) {
  int n = sigma.size();
  arma::mat Sigma = arma::zeros(n, n);
  //TO DO: make this speedier, definitely don't need to compute (i,j) and (j, i) separately
  int w = 0;
  for (unsigned i = 0; i < n; i++) {
    for (unsigned j = 0; j < n; j++) {
      if (i <= j) {
        w = j - i;
      } else {
        w = i - j;
      }
      arma::vec temp = sigma(i)*sigma(j)*exp(-w/kappa);
      Sigma(i, j) = temp(0);
    }
  }
  return(Sigma);
}


//main function
//[[Rcpp::export]]
List RunMCMC(int nIter, int space, arma::cube x, arma::vec delta,
             arma::vec mu, arma::vec mu0, arma::mat Lambda0, arma::mat muPropSigma,
             arma::vec Cx, arma::vec Cx0, arma::vec sigmaX0,
             arma::vec Cy, arma::vec Cy0, arma::vec sigmaY0,
             arma::vec kappa, arma::vec betaKappa0,
             arma::vec sigmaY, arma::vec betaSigmaY0,
             arma::vec thetas) {
  Rcout << "starting " << 0 << std::endl;

  //constants
  int p = mu.size();
  int nStore = nIter/space;

  //functionals
  arma::mat Sigma = compSigma(sigmaY, kappa);
  Rcout << "Finished functionals " << 0 << std::endl;

  //convert data to other forms
  List temp = XToWY(x, Cx, Cy, thetas);
  Rcout << "Finished xToY " << 0 << std::endl;
  
  arma::mat w = temp["w"];
  arma::mat y = temp["y"];

  Rcout << "Finished xToWY " << 0 << std::endl;
  arma::mat SigmaInv = inv(Sigma);
  arma::mat Lambda0Inv = inv(Lambda0);
  Rcout << "Finished data transforms " << 0 << std::endl;


  //storage vectors and matrices
  arma::mat muStore = arma::zeros(p, nStore);
  arma::vec CxStore = arma::zeros(nStore);
  arma::vec CyStore = arma::zeros(nStore);
  arma::vec kappaStore = arma::zeros(nStore);
  arma::mat sigmaYStore = arma::zeros(p, nStore);
  Rcout << "Finished storage " << 0 << std::endl;

  //acceptance rates
  double muRate = 0.0;
  double CxRate = 0.0;
  double CyRate = 0.0;
  double kappaRate = 0.0;
  double sigmaYRate = 0.0;
  Rcout << "Finished acc rates " << 0 << std::endl;


  //quantities used inside the loop (only generate instance once)
  //sampled parameters
  arma::vec muProp(p);
  arma::vec CxProp(1);
  arma::vec CyProp(1);
  arma::vec kappaProp(1);
  arma::vec sigmaProp(p);
  //functional parameters
  arma::mat SigmaProp(p, p);
  //acceptance rate
  double logR(1);
  Rcout << "Finished quants in loop " << 0 << std::endl;


  //Metropolis step loops
  for (unsigned i = 0; i < nIter; i++) {


    ////////////////update mu/////////////////
    muProp = mvrnormCpp(mu, muPropSigma);
    arma::uvec below = arma::find(muProp < 0);
    if (below.size() > 0) {
      logR = -1e16;
    } else {
      logR = (-.5*sumQFCentSq(y, muProp, SigmaInv)
              -.5*sumQFCentSq(muProp, mu0, Lambda0Inv)
              +.5*sumQFCentSq(y, mu, SigmaInv)
              +.5*sumQFCentSq(mu, mu0, Lambda0Inv));
    }
    if (logAccept(logR)) {
      Rcout << "in accept" << i << std::endl;
      mu = muProp;
      muRate++;
    }
    muStore.col(i) = mu;
    Rcout << "stored mu" << i << std::endl;
  }


  //update acceptance rates
  muRate = muRate/nIter;
  CxRate = CxRate/nIter;
  CyRate = CyRate/nIter;
  kappaRate = kappaRate/nIter;
  sigmaYRate = sigmaYRate/nIter;

  //return
  List res;
  res["mu"] = muStore; res["muRate"] = muRate;
  res["Cx"] = CxStore; res["CxRate"] = CxRate;
  res["Cy"] = CyStore; res["CyRate"] = CyRate;
  res["kappa"] = kappaStore; res["kappaRate"] = kappaRate;
  res["sigmaY"] = sigmaYStore; res["sigmaYRate"] = sigmaYRate;
  res["w"] = w; res["y"] = y;

  return(res);
}


