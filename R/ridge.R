#' Solve ridge regression problems
#'
#' @param y The target value - a vector of type numeric.
#' @param X The design matrix - this will not be modified in any way.
#' @param lambda The ridge penalty parameter - set to NULL to try to estimate.
#' @param family Technically allows non-families (such as hurdle/zero-inflated options).
#' Specifies the distribution of your target.
#' @param ... Additional parameters passed to solvers (currently not implemented).
#' @details No choice of link functions or offsets is provided here - all links are canonical,
#' all offsets can be provided within **X** (though whether this makes sense is left to the user to decide -
#' we recommend **not** using this library if you want to use an offset.)
#' @export
ridge <- function( y, X, lambda = 1e-8, family = c("gaussian","poisson","hurdle","zero-inflated"),... ) {

  if( length(family) > 1) {
    family <- family[1]
    rlang::warn( glue::glue("Multiple arguments provided to 'family' - using the first one, '{family}'." ) )
  }

  if( !(family %in% c("gaussian","poisson","hurdle","zero-inflated"))) {
    rlang::abort( glue::glue( "Family {family} not supported." ) )
  }

  solver <- list( gaussian = gaussian_ridge_ls,
                  poisson = poisson_ridge_wls)[[family]]
  do.call(solver, c(list( y = y, X = X, reg_lambda = lambda), list(...) ))
}
