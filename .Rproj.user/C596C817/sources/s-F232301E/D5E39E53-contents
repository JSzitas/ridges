
poisson_bin_hurdle_iwls <- function(y, X, reg_lambda = 1e-8, iter_max = 200, tol = 1e-8,... ) {

  # hurdle - a model on whether something is zero
  is_zero <- binomial_ridge_wls( y = y,
                                 X = X,
                                 reg_lambda = reg_lambda,
                                 iter_max = iter_max,
                                 tol = tol,
                                 ...)

  # counts - only the nonzero parts
  nonzero_indices <- which( y > 0 )
  counts <- poisson_ridge_wls( y = y[nonzero_indices],
                               X = X[nonzero_indices,],
                               reg_lambda = reg_lambda,
                               iter_max = iter_max,
                               tol = tol,
                               ... )
  structure( list( binomial = is_zero,
                   poisson = counts ),
             class = c( "inflated_ridge",
                        "hurdle")
             )
}

poisson_bin_zero_inflated_iwls <- function(y, X, reg_lambda = 1e-8, iter_max = 200, tol = 1e-8,...) {

  # a model on the 'first kind' of zero
  is_zero <- binomial_ridge_wls( y = y,
                                 X = X,
                                 reg_lambda = reg_lambda,
                                 iter_max = iter_max,
                                 tol = tol,
                                 ...)
  # counts - both zero and non-zero parts - this has the 'second kind' of zero
  counts <- poisson_ridge_wls( y = y,
                               X = X,
                               reg_lambda = reg_lambda,
                               iter_max = iter_max,
                               tol = tol,
                               ... )
  structure( list( binomial = is_zero,
                   poisson = counts ),
             class = c("inflated_ridge",
                       "zero_inflated")
             )
}

#' @importFrom stats predict
#' @export
#' @rdname prediction_generics
predict.inflation_model <- function(object, new_data, ...) {

  preds <- purrr::map( object, function(model){
    predict(model, new_data, ...)
  })
  purrr::reduce( preds, `*` )
}

