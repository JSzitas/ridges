# test ridge solver - and optimal lambda estimator

remove(list=ls())
pkgload::load_all()
library(magrittr)

data("CreditCard", package = "AER" )

X <- model.matrix( active~.-1, CreditCard )
y <- CreditCard$active

# set.seed(1071)

rmse <- function( y, pred ) {
  sqrt( mean( (y - pred)^2))
}


lambdas <- c(0,1,-2)
n_folds = 100

glm_rmse <- function( train_y, test_y, train_X, test_X ) {

  df <- data.frame( y = train_y, train_X )
  fit <- glm( y~.-1, data = df, family = poisson )

  pred <- c(predict( fit, as.data.frame(test_X), type = "response"))
  rmse( test_y, pred )
}

rmse_list <- list()

for( i in seq_len(n_folds) ) {
  size <- nrow(X)

  train <- sample( seq_len(size), size * 0.7 )
  test <- seq_len(size)[ -train ]

  train_x <- X[ train, ]
  test_x <- X[ test, ]

  train_y <- y[train]
  test_y <- y[test]

  rmse_list[[i]] <- purrr::map( lambdas, function(lambda) {
    if( lambda < 0 ){
      lambda <- NULL
    }
    model <- poisson_ridge_wls(
      y = train_y,
      X = train_x,
      reg_lambda = lambda,
      tol = 1e-4
    )
    pred <- predict(model, test_x)
    rmse( test_y, pred )
  })

  rmse_list[[i]] <- c( rmse_list[[i]],
                       glm_rmse(train_y, test_y, train_x, test_x) )

}

rmse_list %>%
  do.call(rbind,.) %>%
  data.frame %>%
  lapply( function(i) mean(unlist(i)))

