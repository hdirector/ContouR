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
// [[Rcpp::export]]
double logDet(arma::mat Sigma) {
  return(sum(arma::log(arma::eig_sym(Sigma))));
}



//Compute Sigma based on sigma and kappa
//[[Rcpp::export]]
arma::mat compSigma(arma::vec sigma, arma::vec kappa,  arma::mat thetaDist) {
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





//' main function
// [[Rcpp::export]]
List RunMCMC(int nIter, arma::mat y,
             arma::vec mu, arma::vec mu0, arma::mat Lambda0, arma::vec muPropSD,
             arma::vec kappa,  double betaKappa0, arma::vec kappaPropSD,
             arma::vec sigma, arma::vec betaSigma0, arma::vec sigmaPropSD,
             arma::mat thetaDist) {

  //constants
  int n = y.n_cols;
  int p = mu.size();

  //computed matrices
  arma::mat Sigma = compSigma(sigma, kappa, thetaDist);
  arma::mat SigmaInv = inv(Sigma);
  arma::mat Lambda0Inv = inv(Lambda0);

  //storage vectors and matrices
  arma::mat muStore = arma::zeros(p, nIter);
  arma::vec kappaStore = arma::zeros(nIter);
  arma::mat sigmaStore = arma::zeros(p, nIter);

  //acceptance rates
  arma::vec muRate = arma::zeros(p);
  double kappaRate = 0.0;
  arma::vec sigmaRate = arma::zeros(p);

  //functional parameters
  arma::mat SigmaProp(p, p);
  arma::mat SigmaInvProp(p, p);
  double logR(1);

  //Metropolis step loop
  for (unsigned i = 0; i < nIter; i++) {

    ////////////////update mu/////////////////
    for (int g = 0; g < p; g++) {
      arma::vec muGProp = muPropSD*arma::randn(1) + mu(g);
      arma::vec muProp = mu;
      muProp(g) = muGProp(0);
      logR = (-.5*sumQFCentSq(y, muProp, SigmaInv)
              -.5*sumQFCentSq(muProp, mu0, Lambda0Inv)
              +.5*sumQFCentSq(y, mu, SigmaInv)
              +.5*sumQFCentSq(mu, mu0, Lambda0Inv));
      if (logAccept(logR)) {
        mu = muProp;
        muRate(g) = muRate(g) + 1;
      }
    }
    muStore.col(i) = mu;

    ////////////////update kappa/////////////////
    arma::vec kappaProp =  kappaPropSD*arma::randn(1) + kappa(0);
    SigmaProp = compSigma(sigma, kappaProp, thetaDist);
    SigmaInvProp = inv(SigmaProp);
    if ((kappaProp(0) >= 0) & (kappaProp(0) <= betaKappa0)) {
      logR = (-(n/2)*logDet(SigmaProp) -.5*sumQFCentSq(y, mu, SigmaInvProp)
              +(n/2)*logDet(Sigma) +.5*sumQFCentSq(y, mu, SigmaInv));
      if (logAccept(logR)) {
        kappa = kappaProp;
        kappaRate++;
        Sigma = SigmaProp;
        SigmaInv = SigmaInvProp;
      }
    }
    kappaStore(i) = kappa(0);

    ////////////////update sigmas/////////////////
    for (int g = 0; g < p; g++) {
      arma::vec gProp =  sigmaPropSD*arma::randn(1) + sigma(g);
      if(all(gProp > 0) && all(gProp < betaSigma0(g))) {
        arma::vec sigmaGProp = sigma;
        sigmaGProp(g) = gProp(0);
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
  kappaRate = kappaRate/nIter;
  sigmaRate = sigmaRate/nIter;

  //return
  List res;
  res["mu"] = muStore; res["muRate"] = muRate;
  res["kappa"] = kappaStore; res["kappaRate"] = kappaRate;
  res["sigma"] = sigmaStore; res["sigmaRate"] = sigmaRate;
  return(res);
}
