#include "utils.h"

// [[Rcpp::export]]
List gaussian_ridge( const Eigen::VectorXf& y,
                     const Eigen::MatrixXf& X,
                     double lambda = 1,
                     int iter_max = 100, //this actually does nothing
                     double tolerance = 0.00001) //and neither does this, but they facilitate compatibility
{
  auto tXX = (X.transpose() * X);

  // if negative, find coefficients for lambda == 0 first
  if( lambda < 0 ) {
    lambda = lambda_estimator( y, X, tXX, gaussian_ridge( y, X, 0)['coef'] );
  }

  Eigen::MatrixXf diag_lambda = (Eigen::ArrayXf::Zero(X.cols()) + lambda).matrix().asDiagonal();
  Eigen::VectorXf coef = (tXX + diag_lambda).ldlt().solve(X.transpose() * y);

  return List::create( Named("coef") = coef,
                       Named("lambda") = lambda);
}
