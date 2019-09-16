#ifndef LDSERVER_REGRESSION_H
#define LDSERVER_REGRESSION_H

#include <armadillo>
#include <memory>
#include "Rmath.h"

/**
 * Regression classes ported from rvtests:
 * https://github.com/zhanxw/rvtests
 *
 * Logistic regression original implementation is below. It is updated here to use armadillo and some
 * trivial C++11 features.
 * https://github.com/zhanxw/rvtests/blob/fad75000daf2179f085fa1efc3c86b9a11561d21/regression/LogisticRegression.cpp#L279-L339
 *
 * An earlier non-vectorized version of the logistic regression fitting code found is found below, which is what
 * rvtests was based upon.
 * https://github.com/statgen/bamUtil/blob/c2dd61fe642bca209b0f5be032f50481b990552b/src/LogisticRegression.cpp
 *
 * For background on how the method works, see the following sections in the book "Categorical Data Analysis" by
 * Alan Agresti:
 *   - 4.6: Fitting Generalized Linear Models
 *     - 4.6.1 Newton-Raphson Method
 *   - 5.5: Fitting Logistic Regression Models
 *     - 5.5.4: Newton-Raphson Method Applied to Logistic Regression
 */

class Regression {
public:
  virtual void fit(const arma::vec& y, const arma::mat& x) = 0;
  virtual std::shared_ptr<arma::vec> getResiduals() const = 0;
  virtual std::shared_ptr<arma::vec> getBetas() const = 0;
  virtual std::shared_ptr<arma::vec> getStandardErrors() const = 0;
  virtual std::shared_ptr<arma::mat> getCovBetas() const = 0;
  virtual std::shared_ptr<arma::vec> getPvalues() = 0;
};

class LinearRegression : public Regression {
private:
  arma::vec beta;       // model coefficients, one per column in X
  arma::mat cov_beta;   // asymptotic covariance matrix of betas
  arma::vec pvalue;     // asymptotic wald p-values for each beta

  arma::mat X;          // matrix of covariates, one per column
  arma::mat X_T;        // store transpose of X to avoid re-calculation
  arma::vec Y;          // dependent variable Y
  arma::vec resid;      // residual vector (Y - p) where p is vec of fitted probs
  double sigma2;
  void reset(const arma::mat& x, const arma::vec& y);
public:
  /**
   * Fit linear regression model with OLS.
   * @param y Vector of outcome variable.
   * @param x Matrix of covariates.
   * @param niter Max number of iterations to attempt for convergence.
   */
  void fit(const arma::vec& y, const arma::mat& x) override;
  std::shared_ptr<arma::vec> getPvalues() override;
  std::shared_ptr<arma::vec> getResiduals() const override;
  std::shared_ptr<arma::vec> getBetas() const override;
  std::shared_ptr<arma::vec> getStandardErrors() const override;
  std::shared_ptr<arma::mat> getCovBetas() const override;
  double getSigmaSquared() const;
};

class LogisticRegression : public Regression {
private:
  arma::vec beta;       // model coefficients, one per column in X
  arma::mat cov_beta;   // asymptotic covariance matrix of betas
  arma::vec pvalue;     // asymptotic wald p-values for each beta
  arma::vec p;          // fitted probabilities per y_i
  arma::vec W;          // information matrix I = X' * W * X

  arma::mat X;          // matrix of covariates, one per column
  arma::mat X_T;        // store transpose of X to avoid re-calculation
  arma::vec Y;          // dependent variable Y
  arma::vec eta;        // linear predictor X * beta
  arma::mat info;       // information matrix
  arma::vec delta_beta; // amount to update beta on each iteration
  arma::vec resid;      // residual vector (Y - p) where p is vec of fitted probs
  int niter;            // max number of newton-raphson iterations (default is 100)
  void reset(const arma::mat& x, const arma::vec& y);
public:
  LogisticRegression();
  LogisticRegression(int niter);

  /**
   * Fit a logistic regression model using Newton-Raphson method.
   * @param y Vector of outcome variable.
   * @param x Matrix of covariates.
   * @param niter Max number of iterations to attempt for convergence.
   */
  void fit(const arma::vec& y, const arma::mat& x) override;
  std::shared_ptr<arma::vec> getPvalues() override;
  std::shared_ptr<arma::vec> getPredictedProb() const;
  std::shared_ptr<arma::vec> getResiduals() const override;
  std::shared_ptr<arma::vec> getBetas() const override;
  std::shared_ptr<arma::vec> getStandardErrors() const override;
  std::shared_ptr<arma::mat> getCovBetas() const override;
  double getDeviance();
};

#endif //LDSERVER_REGRESSION_H
