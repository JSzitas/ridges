# test ridge solver - and optimal lambda estimator

remove(list=ls())
pkgload::load_all()
library(magrittr)


df <- rpart::kyphosis

X <- model.matrix( Kyphosis~.-1, df )
y <- as.integer(df$Kyphosis == "present")

set.seed(1071)

acc <- function( y, pred, threshold = 0.5 ) {
  tbl <- table( y, as.integer(pred > threshold))
  sum(diag(tbl))/sum(tbl)
}


lambdas <- c(0,1,-2)
n_folds = 100

glm_acc <- function( train_y, test_y, train_X, test_X ) {

  df <- data.frame( y = train_y, train_X )
  fit <- glm( y~.-1, data = df, family = binomial )

  pred <- c(predict( fit, as.data.frame(test_X), type = "response"))
  acc( test_y, pred )
}

metric_list <- list()

for( i in seq_len(n_folds) ) {
  size <- nrow(X)

  train <- sample( seq_len(size), size * 0.7 )
  test <- seq_len(size)[ -train ]

  train_x <- X[ train, ]
  test_x <- X[ test, ]

  train_y <- y[train]
  test_y <- y[test]

  metric_list[[i]] <- purrr::map( lambdas, function(lambda) {
    if( lambda < 0 ){
      lambda <- NULL
    }
    model <- binomial_ridge_wls(
      y = train_y,
      X = train_x,
      reg_lambda = lambda,
      tol = 1e-4
    )
    pred <- predict(model, test_x)
    acc( test_y, pred )
  })

  metric_list[[i]] <- c( metric_list[[i]],
                         glm_acc(train_y, test_y, train_x, test_x) )

}

metric_list %>%
  do.call(rbind,.) %>%
  data.frame %>%
  lapply( function(i) mean(unlist(i)))
