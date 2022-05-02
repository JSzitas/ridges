poisson_ridge_wls <- function(y, X, reg_lambda = 1e-8, iter_max = 200, tol = 1e-9) {

  n <- nrow(X)
  p <- ncol(X)

  # if null, find coefficients for lambda == 0 first
  if( is.null(reg_lambda) ) {
    ml_coef <- poisson_ridge_wls( y, X, reg_lambda = 0, iter_max, tol )
    tXWX <- ml_coef[["tXWX"]]
    ml_coef <- ml_coef[["coef"]]
    tXWX_eigen <- eigen( tXWX, symmetric = TRUE )
    lambdas <- tXWX_eigen[["values"]]
    # find an estimate of alpha and multiply by coefficients found for unpenalized regression
    alpha <- ml_coef %*% tXWX_eigen[["vectors"]]
    reg_lambda <- stats::median(lambdas/((n-p) + (lambdas* (alpha^2))))
  }

  u <- y + stats::runif(length(y))
  W <- diag(1, nrow = nrow(X))
  z <- log(u) + ((y - u) / u)
  # reg_lambda <- diag(reg_lambda, nrow = ncol(X))

  coef <- rep(0, ncol(X))

  loss <- 0
  iter <- 1

  while (TRUE) {
    if (iter > iter_max) {
      break
    }
    # compute ths ahead of time
    tXWX <- t(X) %*% W %*% X

    new_coef <- solve(tXWX + diag(reg_lambda, nrow = ncol(X))) %*% t(X) %*% W %*% z
    u <- c(exp(X %*% new_coef))
    z <- log(u) + ((y - u) / u)
    W <- diag( u )

    loss_new <- sum(u - y * log(u))
    if (abs(loss_new - loss) < tol) {
      break
    }
    loss <- loss_new
    coef <- new_coef
    iter <- iter + 1
  }

  result <- list(coef = c(coef),
                 lambda = reg_lambda,
                 tXWX = tXWX)

  return(structure( result,
                    class = "poisson_ridge"
  )
  )
}
#' @importFrom stats predict
#' @export
#' @rdname prediction_generics
predict.poisson_ridge <- function(object, new_data, ...) {
  c(exp(new_data %*% object[["coef"]]))
}
