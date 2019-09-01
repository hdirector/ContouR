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
List XToWY(arma::cube x, arma::vec Cx, arma::vec Cy, arma::vec theta) {
  int p = x.n_cols;
  int nSamp = x.n_slices;
  arma::mat w(p, nSamp);
  arma::mat y(p, nSamp);
  arma::vec temp(2);
  arma::vec ptCurr(2); //TO DO, write this in a less dumb way
  arma::vec thetaCurr(1);
  for (unsigned i = 0; i < p; i++) {
    thetaCurr(0) = theta(i);
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
  double pi = 3.14159;
  int n = sigma.size();
  double space = 2*pi/n;
  arma::mat Sigma = arma::zeros(n, n);
  //TO DO: make this speedier, definitely don't need to compute (i,j) and (j, i) separately
  double w = 0;
  for (unsigned i = 0; i < n; i++) {
    for (unsigned j = 0; j < n; j++) {
      if (i <= j) {
        w = std::min(j - i, n - j + i);
      } else {
        w = std::min(i - j, n - i + j);
      }
      w = w*space;
      arma::vec temp = sigma(i)*sigma(j)*exp(-w/kappa);
      Sigma(i, j) = temp(0);
    }
  }
  return(Sigma);
}


//main function
//[[Rcpp::export]]
List RunMCMC(int nIter, arma::cube x, 
             arma::vec mu, arma::vec mu0, arma::mat Lambda0, double muPropSD,
             double nu, 
             arma::vec Cx, double Cx0, double sigmaX0, arma::vec CxPropSD,
             arma::vec Cy, double Cy0, double sigmaY0, arma::vec CyPropSD,
             arma::vec kappa, double alphaKappa0, double betaKappa0, arma::vec kappaPropSD,
             arma::vec sigmaY, double betaSigmaY0, arma::vec sigmaYPropSD,
             arma::vec theta, double theta1PropSD,
             arma::mat kernHat) {
  //Rcout << "starting " << 0 << std::endl;
  
  //constants
  int n = x.n_slices;
  Rcout << "n " << n << std::endl;
  int p = mu.size();
  
  
  //values related to theta
  arma::vec theta1(1);
  theta1(0) = theta(0);
  double thetaSpace = theta(1) - theta(0);
  //Rcout << "a " << 0 << std::endl;
  double theta1UB = theta1(0) + thetaSpace/2;
  //Rcout << "b" << 0 << std::endl;
  

  //functionals
  arma::mat Sigma = compSigma(sigmaY, kappa);
  //Rcout << "Finished functionals " << 0 << std::endl;
  
  //convert data to other forms
  List temp = XToWY(x, Cx, Cy, theta);
  arma::mat w = temp["w"];
  double wSum = arma::accu(w);
  arma::mat y = temp["y"];
  arma::mat SigmaInv = inv(Sigma);
  arma::mat Lambda0Inv = inv(Lambda0);
  arma::vec C(2);
  C(0) = Cx(0);
  C(1) = Cy(0);
  
  //storage vectors and matrices
  arma::mat muStore = arma::zeros(p, nIter);
  arma::mat theta1Store = arma::zeros(nIter);
  arma::vec CxStore = arma::zeros(nIter);
  arma::vec CyStore = arma::zeros(nIter);
  arma::vec kappaStore = arma::zeros(nIter);
  arma::mat sigmaYStore = arma::zeros(p, nIter);
  //Rcout << "Finished storage " << 0 << std::endl;
  
  //acceptance rates
  arma::vec muRate = arma::zeros(p);
  double CxRate = 0.0;
  double CyRate = 0.0;
  double kappaRate = 0.0;
  arma::vec sigmaYRate = arma::zeros(p);
  double theta1Rate = 0.0;
  //Rcout << "Finished acc rates " << 0 << std::endl;
  
  
  //quantities used inside the loop (only generate instance once)
  //sampled parameters
  arma::vec muProp(p);
  arma::vec muPropJ(1);
  arma::vec theta1Prop(1);
  arma::vec thetaProp(p);
  arma::vec thetaPropShift(1);
  arma::mat wPropThetaJ;
  arma::mat yPropThetaJ;
  arma::vec CxProp(1);
  arma::vec CyProp(1);
  arma::vec CProp(2);
  arma::vec kappaProp(1);
  arma::vec sigmaYProp(p);
  arma::vec sigmaYPropJ(1);
  double wSumProp;
  double wSumPropTheta;
  //functional parameters
  arma::mat SigmaProp(p, p);
  arma::mat SigmaInvProp(p, p);
  //acceptance rate
  double logR(1);
  //Rcout << "Finished quants in loop " << 0 << std::endl;
  
  
  //Metropolis step loop
  for (unsigned i = 0; i < nIter; i++) {
    
    ////////////////update mu/////////////////
    for (unsigned j = 0; j < p; j++) {
      muProp = mu;
      muPropJ = muPropSD*arma::randn(1) + mu(j);
      muProp(j) = muPropJ(0);                     
      if(muPropJ(0) < 0) {
        logR = -1e16;
      } else {
        logR = (-.5*sumQFCentSq(y, muProp, SigmaInv)
                -.5*sumQFCentSq(muProp, mu0, Lambda0Inv)
                +.5*sumQFCentSq(y, mu, SigmaInv)
                +.5*sumQFCentSq(mu, mu0, Lambda0Inv));
      }
      if (logAccept(logR)) {
        //Rcout << "in accept" << i << std::endl;
        mu = muProp;
        muRate(j) = muRate(j) + 1;
      }
      muStore(j, i) = mu(j);
    //Rcout << "stored mu" << i << std::endl;
    }
    
    ///////////////update theta///////
    theta1Prop = theta1PropSD*arma::randn(1) + theta1(0);
    thetaPropShift = theta1Prop - theta1;
    thetaProp = theta + thetaPropShift(0);
    temp = XToWY(x, Cx, Cy, thetaProp);
    arma::mat wPropTheta = temp["w"];
    wSumPropTheta = arma::accu(wPropTheta);
    arma::mat yPropTheta = temp["y"];
    if (theta1Prop(0) < 0) {
        logR = -1e16;
    } else if (theta1Prop(0) > theta1UB){
        logR = -1e16;
    } else {
      logR = (-.5*sumQFCentSq(yPropTheta, mu, SigmaInv)
             -(1/(2*pow(nu, 2)))*wSumPropTheta
             +.5*sumQFCentSq(y, mu, SigmaInv)
             +(1/(2*pow(nu, 2)))*wSum);
    }    
    if (logAccept(logR)) {
      theta = thetaProp;
      y = yPropTheta;
      w = wPropTheta;
      wSum = wSumProp;
      theta1 = theta1Prop;
      theta1Rate++;
    }
    theta1Store(i) = theta1(0);


    ////////////////update Cx/////////////////
    CxProp = CxPropSD*arma::randn(1) + Cx;
    CProp = C;
    CProp(0) = CxProp(0);
    temp = XToWY(x, CxProp, Cy, theta);
    arma::mat wPropCx = temp["w"];
    wSumProp = arma::accu(wPropCx);
    arma::mat yPropCx = temp["y"];
    if (!ptInPoly(kernHat, CProp)) { //to do: make ptInPoly take in Cx and Cy seperately
      logR = -1e16;
    } else {
      //Rcout << "in Poly " << i << std::endl;
      logR = (-.5*sumQFCentSq(yPropCx, mu, SigmaInv) 
                -(1/(2*pow(nu, 2)))*wSumProp
                -(1/(2*pow(sigmaX0, 2)))*(CxProp(0) - Cx0)
                +.5*sumQFCentSq(y, mu, SigmaInv)
                +(1/(2*pow(nu, 2)))*wSum
                -(1/(2*pow(sigmaX0, 2)))*(Cx(0) - Cx0));
                //Rcout << "logR " << logR << std::endl;
    }
    if (logAccept(logR)) {
      Cx = CxProp;
      C = CProp;
      w = wPropCx;
      wSum = wSumProp;
      y = yPropCx;
      CxRate++;
    }
    CxStore(i) = Cx(0);
    //Rcout << "stored Cx" << i << std::endl;
    
    ////////////////update Cy/////////////////
    CyProp = CyPropSD*arma::randn(1) + Cy;
    CProp = C;
    CProp(1) = CyProp(0);
    temp = XToWY(x, Cx, CyProp, theta);
    arma::mat wPropCy = temp["w"];
    wSumProp = arma::accu(wPropCy);
    arma::mat yPropCy = temp["y"];
    if (!ptInPoly(kernHat, CProp)) { //to do: make ptInPoly take in Cx and Cy seperately
      logR = -1e16;
    } else {
      logR = (-.5*sumQFCentSq(yPropCy, mu, SigmaInv) 
                -(1/(2*pow(nu, 2)))*wSumProp
                -(1/(2*pow(sigmaY0, 2)))*(CyProp(0) - Cy0)
                +.5*sumQFCentSq(y, mu, SigmaInv)
                +(1/(2*pow(nu, 2)))*wSum
                -(1/(2*pow(sigmaY0, 2)))*(Cy(0) - Cy0));
    }
    if (logAccept(logR)) {
      Cy = CyProp;
      C = CProp;
      w = wPropCy;
      wSum = wSumProp;
      y = yPropCy;
      CyRate++;
    }
    CyStore(i) = Cy(0);
    //Rcout << "stored Cy" << i << std::endl;
    
    
    ////////////////update kappa///////////////// 
    kappaProp =  kappaPropSD*arma::randn(1) + kappa;
    if ((kappaProp(0) < alphaKappa0) | (kappaProp(0) > betaKappa0)) {
      logR = -1e16;
    } else {
      SigmaProp = compSigma(sigmaY, kappaProp);
      SigmaInvProp = inv(SigmaProp);
      logR = (-(n/2)*logDet(SigmaProp) + 
        -.5*sumQFCentSq(y, mu, SigmaInvProp) 
        +(n/2)*logDet(Sigma) + 
        +.5*sumQFCentSq(y, mu, SigmaInv)); 
    }
    if (logAccept(logR)) {
      kappa = kappaProp;
      kappaRate++;
    }
    kappaStore(i) = kappa(0);
    
    ////////////////update sigmas/////////////////
    for (unsigned j = 0; j < p; j++) {
      sigmaYProp = sigmaY;
      sigmaYPropJ = sigmaYPropSD*arma::randn(1) + sigmaY(j);
      if ((sigmaYPropJ(0) < 0) | (sigmaYPropJ(0) > betaSigmaY0)) {
        logR = -1e16;
      } else {
        sigmaYProp(j) = sigmaYPropJ(0);
        SigmaProp = compSigma(sigmaYProp, kappa);
        SigmaInvProp = inv(SigmaProp);
        logR = (-(n/2)*logDet(SigmaProp) +
          -.5*sumQFCentSq(y, mu, SigmaInvProp)
          +(n/2)*logDet(Sigma) +
          +.5*sumQFCentSq(y, mu, SigmaInv));
      }
      if (logAccept(logR)) {
        sigmaY = sigmaYProp;
        Sigma = SigmaProp;
        SigmaInv = SigmaInvProp;
        sigmaYRate(j) = sigmaYRate(j) + 1;
      }
      sigmaYStore(j, i) = sigmaY(j);
    }
  }
  
  
  
  //update acceptance rates
  muRate = muRate/nIter;
  theta1Rate = theta1Rate/nIter;
  CxRate = CxRate/nIter;
  CyRate = CyRate/nIter;
  kappaRate = kappaRate/nIter;
  sigmaYRate = sigmaYRate/nIter;
  
  //return
  List res;
  res["mu"] = muStore; res["muRate"] = muRate;
  res["theta1"] = theta1Store; res["thetaRate"] = theta1Rate;
  res["Cx"] = CxStore; res["CxRate"] = CxRate;
  res["Cy"] = CyStore; res["CyRate"] = CyRate;
  res["kappa"] = kappaStore; res["kappaRate"] = kappaRate;
  res["sigmaY"] = sigmaYStore; res["sigmaYRate"] = sigmaYRate;
  res["w"] = w; res["y"] = y;
  
  return(res);
}


