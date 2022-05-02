#ifndef util_header
#define util_header

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <algorithm>
using namespace Rcpp;
using namespace RcppEigen;
using namespace Eigen;

double median( Eigen::VectorXf x )
{
  std::sort(x.data(),x.data()+x.size());
  return x.size() % 2 == 0 ?
  x.segment( (x.size()-2)/2, 2 ).mean() :
    x( x.size()/2 );
}

// template<typename Function>
double lambda_estimator( const Eigen::VectorXf& y,
                         const Eigen::MatrixXf& X,
                         const Eigen::MatrixXf& tXX,
                         const Eigen::VectorXf ml_coef,
                         // Function ridge_estimator,
                         int iter_max = 100,
                         double tolerance = 0.00001
) {
    auto n = X.rows();
    auto p = X.cols();

    // Eigen::VectorXf ml_coef = ridge_estimator( y, X, 0, iter_max, tolerance )["coef"];
    // the weighed covariance matrix is symmetric
    SelfAdjointEigenSolver<MatrixXf> tXX_eigen(tXX);
    // grab eigenvectors
    auto eigenvecs = tXX_eigen.eigenvectors();
    // grab eigenvalues
    auto eigenvals = tXX_eigen.eigenvalues();
    // find an estimate of alpha and multiply by coefficients found for un-penalized regression
    auto alpha = (ml_coef.transpose() * eigenvecs).transpose();

    auto res_lambdas = eigenvals.array() *
      (1/((n-p) + eigenvals.array() * (alpha.array() * alpha.array()) ));
    double lambda = median( res_lambdas );

  return lambda;
}

#endif
