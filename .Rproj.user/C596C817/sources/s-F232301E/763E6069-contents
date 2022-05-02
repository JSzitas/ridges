
logit_link <- function( x ) {
  y <- c(1/(1+exp(-x)))
  # y[ y < 0.0000001 ] <- 0.0000001
  # y[ y > 0.9999999 ] <- 0.9999999
  y[ y < 0.001 ] <- 0.001
  y[ y > 0.999 ] <- 0.999
  return(y)
}

# logit_link <- function(x) {
#   exp(x)/(1 + exp(x))
# }


binomial_ridge_wls <- function(y, X, reg_lambda = 1e-8, iter_max = 300, tol = 1e-9) {

  estimate_lambda <- FALSE
  # if null, find coefficients for lambda == 0 first
  if( is.null(reg_lambda) ) {
    estimate_lambda <- TRUE
    reg_lambda <- 0
    ml_coef <- binomial_ridge_wls( y, X, reg_lambda = 0, iter_max, tol )[["coef"]]
  }

  u <- logit_link( X %*% rep( 0,ncol(X)))
  W <- diag(1, nrow = nrow(X))
  # return(list( y, u))
  z <- (y - u) / ( u * (1-u))
  reg_lambda <- diag(reg_lambda, nrow = ncol(X))

  coef <- rep(0, ncol(X))

  loss <- 0
  iter <- 1

  n <- nrow(X)
  p <- ncol(X)

  while (TRUE) {
    if (iter > iter_max) {
      break
    }

    tXWX <- t(X) %*% W %*% X

    # this allows us to compute the ridge lambda parameter, using a whack estimator from
    # "Poisson Ridge Regression Estimators: A Performance Test" (2021),
    # known in that paper as "k36" (k37 looks potentially promising but, adding a potentially confusing choice
    # is not something I am too keen on, as it is unclear exactly when which estimator is better - and it is
    # used here primarily to avoid requiring more confusing choices from the user. )
    if( estimate_lambda ) {
      # the weighed covariance matrix is symmetric
      tXWX_eigen <- eigen( tXWX, symmetric = TRUE )
      lambdas <- tXWX_eigen[["values"]]
      # find an estimate of alpha and multiply by coefficients found for unpenalized regression
      alpha <- ml_coef %*% tXWX_eigen[["vectors"]]
      reg_lambda <- stats::median( lambdas/((n-p) + (lambdas* (alpha^2))))
      reg_lambda <- diag(reg_lambda, nrow = ncol(X))
    }

    new_coef <- solve(tXWX + reg_lambda) %*% t(X) %*% z
    u <- c(logit_link( c(X %*% new_coef) ))
    z <- (y - u)
    u <- u * (1-u)
    W <- diag( u )

    # loss_new <- sum( -y*log(u) -(1-y)*log(u) )
    if ( all(abs(new_coef - coef) < tol ) ) {
      break
    }
    # loss <- loss_new
    coef <- new_coef
    iter <- iter + 1
  }
  return(structure( list(coef = c(coef),
                         lambda = reg_lambda),
                    class = "logistic_ridge"
  ))
}
#' #' @importFrom stats predict
#' @export
#' @rdname prediction_generics
predict.logistic_ridge <- function(object, new_data, ...) {
  logit_link(new_data %*% object[["coef"]])
}

