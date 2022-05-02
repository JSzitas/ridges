
# this does not require IWLS :)
gaussian_ridge_ls <- function( y, X, reg_lambda = 1e-8, ... ) {

  tXX <- t(X) %*% X
  # if null, find coefficients for lambda == 0 first
  if( is.null(reg_lambda) ) {
    ml_coef <- gaussian_ridge_ls( y, X, reg_lambda = 0 )[["coef"]]
    # the weighed covariance matrix is symmetric
    tXX_eigen <- eigen( tXX, symmetric = TRUE )
    lambdas <- tXX_eigen[["values"]]
    # find an estimate of alpha and multiply by coefficients found for unpenalized regression
    alpha <- ml_coef %*% tXX_eigen[["vectors"]]
    n <- nrow(X)
    p <- ncol(X)
    reg_lambda <- stats::median( lambdas/((n-p) + (lambdas* (alpha^2))))
  }

  coef <- y %*% X %*% solve(tXX + (reg_lambda * diag(ncol(X))))

  return(structure( list(coef = c(coef),
                         lambda = reg_lambda),
                    class = "gaussian_ridge" )
  )
}
#' #' @importFrom stats predict
#' @export
#' @rdname prediction_generics
#' @param object The fitted model.
#' @param ... Additional parameters (only for generic consistency).
#' @param new_data New data (always a matrix).
#' @return Predictions.
predict.gaussian_ridge <- function( object, new_data, ... ) {
  c(new_data %*% object[["coef"]])
}
