#include "Regression.h"
using namespace arma;
using namespace std;

void LogisticRegression::reset(const mat& x, const vec& y) {
  int nrows = x.n_rows;
  int ncols = x.n_cols;

  beta.zeros(ncols);
  cov_beta.zeros(ncols, ncols);
  pvalue.zeros(ncols);
  p.zeros(nrows);
  W.zeros(nrows);
  info.zeros(nrows, nrows);
  delta_beta.zeros(ncols);

  X = mat(x);
  X_T = arma::trans(X);
  Y = vec(y);
}

double LogisticRegression::getDeviance() {
  return -2.0 * arma::sum((Y % arma::log(p)) + ((1 - Y) % arma::log(1 - p)));
}

void LogisticRegression::fit(const vec& y, const mat& x, int niter) {
  // Check number of iterations.
  if (niter <= 0) {
    throw std::invalid_argument("niter must be > 0");
  }

  // Clear previous estimates.
  reset(x, y);

  // Begin Newton-Raphson iterations.
  int rounds = 0;
  double lastDeviance = -99999;
  double currentDeviance;
  arma::vec XT_resid;
  while (rounds < niter) {
    eta = x * beta;
    p = 1 / (1 + exp(-eta));
    W = p % (1.0 - p);                // The % gets arma to multiply each element w/o transpose
    info = X_T * diagmat(W) * X;      // This is the current estimate of the information matrix
    XT_resid = X_T * (Y - p);
    delta_beta = solve(info, XT_resid);  // This solves (X'WX)^-1 * delta_beta = X'(y - p) for delta_beta
    beta += delta_beta;
    currentDeviance = getDeviance();

    if (rounds > 1 && fabs(currentDeviance - lastDeviance) < 1e-3) {
      rounds = 0;
      break;
    }

    if (fpclassify(currentDeviance) != FP_NORMAL) {
      throw runtime_error("Deviance became subnormal or over/underflowed during LR fit");
    }

    lastDeviance = currentDeviance;
    rounds++;
  }

  if (rounds == niter) {
    throw runtime_error("Not enough iterations during LR fit");
  }

  // Compute asymptotic covariance matrix of betas (inverse of information matrix).
  cov_beta = arma::inv_sympd(info);

  // Compute residuals (these are "response" residuals, not deviance residuals).
  resid = Y - p;
}

shared_ptr<vec> LogisticRegression::getPvalues() {
  if (beta.n_elem <= 0 || sum(beta) == 0.0) {
    throw runtime_error("Call fit() before getting p-values for logistic regression");
  }

  for (int j = 0; j < beta.n_elem; j++) {
    double zstat = beta[j] * beta[j] / cov_beta(j, j);
    pvalue[j] = pchisq(zstat, 1.0, 0, 0);
  }

  auto ret = make_shared<vec>(pvalue);
  return ret;
}

std::shared_ptr<arma::vec> LogisticRegression::getPredictedProb() const {
  return make_shared<vec>(p);
}

std::shared_ptr<arma::vec> LogisticRegression::getResiduals() const {
  return make_shared<vec>(resid);
}

std::shared_ptr<arma::vec> LogisticRegression::getBetas() const {
  return make_shared<vec>(beta);
}

std::shared_ptr<arma::vec> LogisticRegression::getStandardErrors() const {
  return make_shared<vec>(arma::sqrt(diagvec(cov_beta)));
}