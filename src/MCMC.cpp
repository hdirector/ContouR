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

// compute log of the determinant using that the product of
// the eigenvalues of a matrix equals its determinant
double logDet(arma::mat Sigma) {
  return(sum(arma::log(arma::eig_sym(Sigma))));
}


// find the distance between two points
// [[Rcpp::export]]
double ptDist(arma::vec p1, arma::vec p2) {
  return(sqrt(arma::sum(arma::square(p1 - p2))));
}

//check if a point is in a polygon using the ray-casting approach
//e.g., https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
// [[Rcpp::export]]
bool ptInPoly(arma::mat poly, double ptX, double ptY) {
  int n = poly.n_cols;
  int count = 0; //count of intersection with test ray and polygon edges
  
  //main section: loop over the edges and count intersections with test ray
  
  
  for (unsigned i = 0; i < n - 1; i++) {
    arma::vec ea = poly.col(i);
    arma::vec eb = poly.col(i + 1);
    //horizontal edge
    if (ea(1) == eb(1)) {
      if (ea(1) == ptY) { //test point and edge must be at same height for point to be in
        if (ptX >= std::min(ea(0), eb(0)) & ptX <= std::max(ea(0), eb(0))) {
          return(TRUE);
        }
      }
      
      //vertical edge
    } else if (ea(0) == eb(0)) {
      if (ea(0) > ptX) { //line must be to the right
        if ((ptY >= std::min(ea(1), eb(1))) &&
            (ptY <= std::max(ea(1), eb(1)))) {//test point y must be between edge end points y's
          count += 1;
        }
      } else if (ea(0) == ptX) {
        if ((ptY >= std::min(ea(1), eb(1))) &&
            (ptY <= std::max(ea(1), eb(1)))) { //test point on boundary line
          return(TRUE);
        }
      }
      
      //typical edge
    } else {
      //find x-value at which the edge intersects the height of the test
      //ray. Then check if that x-value is on or to the right of the point
      double m = (eb(1) - ea(1))/(eb(0) - ea(0));
      double b = eb(1) - m*eb(0);
      double xInter = (ptY - b)/m;
      if (abs(xInter - ptX) <= 1e-10) { //test point is (approximately) on edge
        return(TRUE);
      } else if ((xInter > ptX) && (xInter >= std::min(ea(0), eb(0)))
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
arma::vec projDist(arma::vec Cx, arma::vec Cy, arma::vec theta,
                   double ptX, double ptY) {
  double pi = 3.14159;
  double eps = .00001;
  arma::vec wSq(1);
  arma::vec y(1);
  //vertical line
  if  (((theta(0) >=  -eps) & (theta(0) <=  eps)) |
       ((theta(0) >= pi - eps) & (theta(0) <= pi + eps))) {
    wSq = pow(ptY - Cy, 2);
    y = abs(ptX - Cx);
    //horizontal line
  } else if (((theta(0) >=  pi/2 - eps) & (theta(0) <=  pi/2 + eps)) |
    ((theta(0) >= 3*pi/2 - eps) & (theta(0) <= 3*pi/2 + eps))) {
    wSq = pow(ptX - Cx, 2);
    y = abs(ptY - Cy);
    //typical case
  } else {
    //equation of generating line (line that extends from (Cx, Cy) along theta)
    arma::vec m = sin(theta)/cos(theta);
    arma::vec b = Cy - m*Cx;
    //equation of perpendicular line
    arma::vec mProj = -(1/m);
    arma::vec bProj = ptY - mProj*ptX;
    //intersection of generating line and perpendicular line
    arma::vec xInter = (b - bProj)/(mProj - m);
    arma::vec yInter = m*xInter + b;
    //w value (perpendicular length)
    wSq = arma::square(ptX - xInter) + arma::square(ptY - yInter);
    y = sqrt(arma::square(Cx(0) - xInter) + arma::square(Cy(0) - yInter));
  }

  arma::vec ret(2);
  ret(0) = wSq(0);
  ret(1) = y(0);
  return(ret);
}

//Compute W and Y values
//[[Rcpp::export]]
List XToWY(arma::cube x, arma::vec Cx, arma::vec Cy, arma::vec theta) {
  int p = x.n_cols;
  int nSamp = x.n_slices;
  arma::mat wSq(p, nSamp);
  arma::mat y(p, nSamp);
  arma::vec temp(2);
  arma::vec thetaCurr(1);
  for (unsigned i = 0; i < p; i++) {
    thetaCurr(0) = theta(i);
    for (unsigned j = 0; j < nSamp; j ++) {
      temp = projDist(Cx, Cy, thetaCurr,  x(0,i, j), x(1,i, j));
      wSq(i, j) = temp(0);
      y(i, j) = temp(1);
    }
  }

  List wy;
  wy["wSq"] = wSq;
  wy["y"] = y;
  return(wy);
}

//compute distances between angles
//[[Rcpp::export]]
arma::mat compThetaDist(int n, double space) {
  arma::mat thetaDist = arma::zeros(n, n);
  double w = 0;
  for (unsigned i = 0; i < n; i++) {
    for (unsigned j = 0; j <= i; j++) {
      if (i <= j) {
        w = std::min(j - i, n - j + i);
      } else {
        w = std::min(i - j, n - i + j);
      }
      thetaDist(i, j) = w*space;
      thetaDist(j, i) = w*space;
    }
  }
  return(thetaDist);
}


//Compute Sigma based on sigma and kappa
//[[Rcpp::export]]
arma::mat compSigma(arma::vec sigma, arma::vec kappa, arma::mat thetaDist) {
  int n = sigma.size();
  arma::mat Sigma = arma::zeros(n, n);
  for (unsigned i = 0; i < n; i++) {
    for (unsigned j = 0; j <= i; j++) {
      arma::vec temp = sigma(i)*sigma(j)*exp(-thetaDist(i, j)/kappa);
      Sigma(i, j) = temp(0);
      Sigma(j, i) = temp(0);
    }
  }
  return(Sigma);
}

//main function
//[[Rcpp::export]]
List RunMCMC(int nIter, arma::cube x,
             arma::vec mu, arma::vec mu0, arma::mat Lambda0, arma::mat muPropCov,
             arma::vec nu, double nuPropSD, double alphaNu0, double betaNu0,
             arma::vec Cx, double Cx0, double sigmaX0, arma::vec CxPropSD,
             arma::vec Cy, double Cy0, double sigmaY0, arma::vec CyPropSD,
             arma::vec kappa, double alphaKappa0, double betaKappa0,
             arma::vec kappaPropSD,
             arma::vec sigma, double betaSigma0, arma::mat sigmaPropCov,
             arma::vec theta1, double theta1PropSD,
             arma::mat kernHat,
             arma::uvec gStart, arma::uvec gEnd) {

  //constants
  double pi = 3.14159;
  int n = x.n_slices;
  int p = mu.size();
  int nGroup = gStart.size();

  //compute theta info
  double thetaSpace = 2*pi/p;
  arma::vec seq = arma::linspace(0, p - 1, p);
  arma::vec theta = theta1(0) + seq*thetaSpace;
  double theta1UB = theta1(0) + thetaSpace/2;

  //Compute sigma matrix
  arma::mat thetaDist = compThetaDist(p, thetaSpace);
  arma::mat Sigma = compSigma(sigma, kappa, thetaDist);
  arma::mat SigmaInv = inv(Sigma);

  //convert data to other forms
  List temp = XToWY(x, Cx, Cy, theta);
  arma::mat wSq = temp["wSq"];
  double wSqSum = arma::accu(wSq);
  arma::mat y = temp["y"];
  arma::mat Lambda0Inv = inv(Lambda0);

  //storage vectors and matrices
  arma::mat muStore = arma::zeros(p, nIter);
  arma::mat theta1Store = arma::zeros(nIter);
  arma::vec nuStore = arma::zeros(nIter);
  arma::vec CxStore = arma::zeros(nIter);
  arma::vec CyStore = arma::zeros(nIter);
  arma::vec kappaStore = arma::zeros(nIter);
  arma::mat sigmaStore = arma::zeros(p, nIter);

  //acceptance rates
  arma::vec muRate = arma::zeros(nGroup);
  double theta1Rate = 0.0;
  double nuRate = 0.0;
  double CxRate = 0.0;
  double CyRate = 0.0;
  double kappaRate = 0.0;
  arma::vec sigmaRate = arma::zeros(nGroup);

  //functional parameters
  arma::mat SigmaProp(p, p);
  arma::mat SigmaInvProp(p, p);
  double logR(1);

  //Metropolis step loop
  for (unsigned i = 0; i < nIter; i++) {
    
    ////////////////update mu/////////////////
    for (int g = 0; g < nGroup; g++) {
      arma::uvec gInd = arma::linspace<arma::uvec>(gStart(g), gEnd(g),
                                                   gEnd(g) - gStart(g) + 1);
      arma::vec muGProp = mvrnormCpp(mu(gInd), muPropCov(gInd, gInd));      
      if (all(muGProp > 0)) {
        arma::vec muProp = mu;
        muProp(gInd) = muGProp;
        logR = (-.5*sumQFCentSq(y, muProp, SigmaInv)
                -.5*sumQFCentSq(muProp, mu0, Lambda0Inv)
                +.5*sumQFCentSq(y, mu, SigmaInv)
                +.5*sumQFCentSq(mu, mu0, Lambda0Inv));
        if (logAccept(logR)) {
          mu = muProp;
          muRate(g) = muRate(g) + 1;
        }
      }
    }
    muStore.col(i) = mu;
    

    ///////////////update theta 1///////
    arma::vec theta1Prop = theta1PropSD*arma::randn(1) + theta1(0);
    arma::vec thetaPropShift = theta1Prop - theta1;
    arma::vec thetaProp = theta + thetaPropShift(0);
    temp = XToWY(x, Cx, Cy, thetaProp);
    arma::mat wSqPropTheta = temp["wSq"];
    double wSqSumPropTheta = arma::accu(wSqPropTheta);
    arma::mat yPropTheta = temp["y"];
    if ((theta1Prop(0) >= 0) & (theta1Prop(0) <= theta1UB)) {
      logR = (-.5*sumQFCentSq(yPropTheta, mu, SigmaInv)
             -(1/(2*pow(nu(0), 2)))*wSqSumPropTheta
             +.5*sumQFCentSq(y, mu, SigmaInv)
             +(1/(2*pow(nu(0), 2)))*wSqSum);
      if (logAccept(logR)) {
        theta = thetaProp;
        y = yPropTheta;
        wSq = wSqPropTheta;
        wSqSum = wSqSumPropTheta;
        theta1 = theta1Prop;
        theta1Rate++;
      }
    }
    theta1Store(i) = theta1(0);
    
    ////////////update nu////////////////
    arma::vec nuProp = nuPropSD*arma::randn(1) + nu(0);
    if ((nuProp(0) > alphaNu0) & (nuProp(0) < betaNu0)) {
      logR = (-n*p*log(nuProp(0))/2 -(1/(2*pow(nuProp(0), 2)))*wSqSum
              +n*p*log(nu(0))/2 + (1/(2*pow(nu(0), 2)))*wSqSum);
      if (logAccept(logR)) {
        nu = nuProp;
        nuRate++;
      }
    }
    nuStore(i) = nu(0);


    ////////////////update Cx/////////////////
    arma::vec CxProp = CxPropSD*arma::randn(1) + Cx;
    if (ptInPoly(kernHat, CxProp(0), Cy(0))) {
      temp = XToWY(x, CxProp, Cy, theta);
      arma::mat wSqPropCx = temp["wSq"];
      double wSqSumPropCx = arma::accu(wSqPropCx);
      arma::mat yPropCx = temp["y"];
      logR = (-.5*sumQFCentSq(yPropCx, mu, SigmaInv)
              -(1/(2*pow(nu(0), 2)))*wSqSumPropCx
              -(1/(2*pow(sigmaX0, 2)))*(CxProp(0) - Cx0)
              +.5*sumQFCentSq(y, mu, SigmaInv)
              +(1/(2*pow(nu(0), 2)))*wSqSum
              -(1/(2*pow(sigmaX0, 2)))*(Cx(0) - Cx0));

      if (logAccept(logR)) {
        Cx = CxProp;
        wSq = wSqPropCx;
        wSqSum = wSqSumPropCx;
        y = yPropCx;
        CxRate++;
      }
    }
    CxStore(i) = Cx(0);

    ////////////////update Cy/////////////////
    arma::vec CyProp = CyPropSD*arma::randn(1) + Cy;
    if (ptInPoly(kernHat, Cx(0), CyProp(0))) {
      temp = XToWY(x, Cx, CyProp, theta);
      arma::mat wSqPropCy = temp["wSq"];
      double wSqSumPropCy = arma::accu(wSqPropCy);
      arma::mat yPropCy = temp["y"];
      logR = (-.5*sumQFCentSq(yPropCy, mu, SigmaInv)
              -(1/(2*pow(nu(0), 2)))*wSqSumPropCy
              -(1/(2*pow(sigmaY0, 2)))*(CyProp(0) - Cy0)
              +.5*sumQFCentSq(y, mu, SigmaInv)
              +(1/(2*pow(nu(0), 2)))*wSqSum
              -(1/(2*pow(sigmaY0, 2)))*(Cy(0) - Cy0));
      if (logAccept(logR)) {
        Cy = CyProp;
        wSq = wSqPropCy;
        wSqSum = wSqSumPropCy;
        y = yPropCy;
        CyRate++;
      }
    }
    CyStore(i) = Cy(0);

    ////////////////update kappa/////////////////
    arma::vec kappaProp =  kappaPropSD*arma::randn(1) + kappa;
    if ((kappaProp(0) >= alphaKappa0) & (kappaProp(0) <= betaKappa0)) {
      SigmaProp = compSigma(sigma, kappaProp, thetaDist);
      SigmaInvProp = inv(SigmaProp);
      logR = (-(n/2)*logDet(SigmaProp)  -.5*sumQFCentSq(y, mu, SigmaInvProp)
              +(n/2)*logDet(Sigma)  +.5*sumQFCentSq(y, mu, SigmaInv));
      if (logAccept(logR)) {
        kappa = kappaProp;
        kappaRate++;
      }
    }
    kappaStore(i) = kappa(0);

    ////////////////update sigmas/////////////////
    for (int g = 0; g < nGroup; g++) {
      arma::uvec gInd = arma::linspace<arma::uvec>(gStart(g), gEnd(g),
                                                   gEnd(g) - gStart(g) + 1);
      arma::vec gProp = mvrnormCpp(sigma(gInd), sigmaPropCov(gInd, gInd));
      if (all(gProp >= 0) && all(gProp < betaSigma0)) {
        arma::vec sigmaGProp = sigma; 
        sigmaGProp(gInd) = gProp;
        arma::mat SigmaGProp = compSigma(sigmaGProp, kappa, thetaDist);
        arma::mat SigmaGInvProp = inv(SigmaGProp);
        logR = (-(n/2)*logDet(SigmaGProp) -.5*sumQFCentSq(y, mu, SigmaGInvProp)
                +(n/2)*logDet(Sigma) +.5*sumQFCentSq(y, mu, SigmaInv));
        if (logAccept(logR)) {
          sigma = sigmaGProp;
          Sigma = SigmaGProp;
          SigmaInv = SigmaGInvProp;
          sigmaRate(g) = sigmaRate(g) + 1;
        }
      }
    }
    sigmaStore.col(i) = sigma;
  }
  

  //update acceptance rates
  muRate = muRate/nIter;
  theta1Rate = theta1Rate/nIter;
  nuRate = nuRate/nIter;
  CxRate = CxRate/nIter;
  CyRate = CyRate/nIter;
  kappaRate = kappaRate/nIter;
  sigmaRate = sigmaRate/nIter;

  //return
  List res;
  res["mu"] = muStore; res["muRate"] = muRate;
  res["theta1"] = theta1Store; res["thetaRate"] = theta1Rate;
  res["nu"] = nuStore; res["nuRate"] = nuRate;
  res["Cx"] = CxStore; res["CxRate"] = CxRate;
  res["Cy"] = CyStore; res["CyRate"] = CyRate;
  res["kappa"] = kappaStore; res["kappaRate"] = kappaRate;
  res["sigma"] = sigmaStore; res["sigmaRate"] = sigmaRate;
  res["wSq"] = wSq; res["y"] = y;

  return(res);
}


