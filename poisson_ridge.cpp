#include "utils.h"

// [[Rcpp::export]]
List poisson_ridge( const Eigen::VectorXf& y,
                    const Eigen::MatrixXf& X,
                    double lambda = 1,
                    int iter_max = 100,
                    double tolerance = 0.00001 )
{

  bool estimate_lambda = false;
  if( lambda < 0 ) {
    estimate_lambda = true;
  }

  Eigen::VectorXf u = y;
  Eigen::MatrixXf W = (Eigen::ArrayXf::Zero(X.cols()) + 1.01).matrix().asDiagonal();
  Eigen::VectorXf z = log(u) + ((y - u) / u);

  auto tXWX = X.transpose() * W * X;

  Eigen::VectorXf coef(X.cols());
  Eigen::VectorXf new_coef(X.cols());
  double loss = 0;
  double loss_new = 0;
  int iter = 1;

  while (true) {
    if (iter > iter_max) {
      break
    }
    // compute current lambda
    if( estimate_lambda ) {
      lambda = lambda_estimator( y, X, tXWX, poisson_ridge( y, X, 0)['coef'] );
    }

    new_coef = solve(tXWX + reg_lambda) * X.transpose() * W * z;
    u = exp(X * new_coef);
    z = log(u) + ((y - u) / u);
    W = diag( u );

    loss_new = sum(u - y * log(u));
    if (abs(loss_new - loss) < tol) {
      break
    }
    loss = loss_new;
    coef = new_coef;
    iter = iter + 1;
    tXWX = X.transpose() * W * X;
  }


  return List::create( Named("coef") = coef,
                       Named("lambda") = lambda);
}
