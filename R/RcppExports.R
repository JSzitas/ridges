# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

gaussian_ridge <- function(y, X, lambda = 1, iter_max = 100L, tolerance = 0.00001) {
    .Call('_ridges_gaussian_ridge', PACKAGE = 'ridges', y, X, lambda, iter_max, tolerance)
}

