
#' Bayesian Linear Regression with Normal-Inverse Gamma prior.
#'
#' @param formula A string fromula.
#' @param data A data frame with your data including response.
#' @param m A matrix with the hyperparameter of the mean for the normal distribution.
#' @param M A matrix with the hyperparameter of the covariance for the normal distribution.
#' @param a A number of hyperparameter for the gamma distribution.
#' @param b A number of hyperparameter for the gamma distribution.
#' @return A list with the estimated coefficients and standard deviations of the model from the Classical and Bayesian approach
#' @examples
#' # Noninformative Prior (same results Bayesian and Classical Linear Regression)
#' m <- matrix(0, nrow = 2, ncol = 1)
#' M <- diag(10000000000000000000, 2)
#' a <- 0.000000000000000000000000001
#' b <- 0.000000000000000000000000001
#' test <- posterior_nig(
#' formula = 'dist ~ speed', 
#' data = cars,
#' m = m,
#' M = M,
#' a = a,
#' b = b 
#' )
#' test$classic_coef
#' test$classic_sigma
#' test$bayesian_coef
#' test$bayesian_sigma

posterior_nig <- function(formula, data, m, M, a, b){
  lm_classic <- lm(formula, data)
  X <- model.matrix(lm_classic)
  y <- as.matrix(lm_classic$model[1])
  n <- dim(y)[1]
  M_tilde <- solve(t(X)%*%X + solve(M))
  m_tilde <- M_tilde %*% (solve(M)%*%m + t(X)%*%y)
  a_tilde <- a + n/2
  b_tilde <- b + (1/2)*(
    t(y) %*% y 
    + t(m) %*% solve(M) %*% m 
    - t(m_tilde) %*% solve(M_tilde) %*% m_tilde
    )
  a_tilde <- as.numeric(a_tilde)
  b_tilde <- as.numeric(b_tilde)
  sigma <- sqrt(b_tilde/(a_tilde-1))
  print('Done!')
  mylist <- list(
    classic_coef = lm_classic$coef,
    classic_sigma = summary(lm_classic)$sigma,
    bayesian_coef = m_tilde,
    bayesian_sigma = sigma
  )
  return(mylist)
}

#' Bayesian Linear Regression with  Zellner's prior.
#'
#' @param formula A string fromula.
#' @param data A data frame with your data including response.
#' @param m1 A number of hyperparameter of the intercept.
#' @param g A string of hyperparameter of the Zellner's prior. "typeI" is 1/n, "typeII" is 1/k^2 and "typeIII" is 1/max(n,k^2). 
#' @param a A number of hyperparameter for the gamma distribution.
#' @param b A number of hyperparameter for the gamma distribution.
#' @return A list with the estimated coefficients and standard deviations of the model from the Classical and Bayesian approach
#' @examples
#' a <- 0.0001
#' b <- 0.0001
#' test <- zellner(
#' formula = 'dist ~ speed', 
#' data = cars,
#' m1 = 0,
#' g = 'typeI',
#' a = a,
#' b = b 
#' )
#' test$classic_coef
#' test$classic_sigma
#' test$bayesian_coef
#' test$bayesian_sigma

zellner <- function(formula, data, m1 = m1, g = g, a = a, b = b){
  lm_classic <- lm(formula, data)
  X <- model.matrix(lm_classic)
  y <- as.matrix(lm_classic$model[1])
  n <- dim(y)[1]
  p <- dim(X)[2]
  if(g == 'typeI'){g=1/n}
  else if (g == 'typeII'){g=1/p}
  else if (g == 'typeIII'){g=1/max(c(n,p))}
  M <- solve(g * t(X)%*% X)
  m <- c(m1, rep(0, (p-1)))
  M_tilde <- solve(t(X)%*%X + solve(M))
  m_tilde <- M_tilde %*% (solve(M)%*%m + t(X)%*%y)
  a_tilde <- a + n/2
  b_tilde <- b + (1/2)*(
    t(y) %*% y 
    + t(m) %*% solve(M) %*% m 
    - t(m_tilde) %*% solve(M_tilde) %*% m_tilde
  )
  a_tilde <- as.numeric(a_tilde)
  b_tilde <- as.numeric(b_tilde)
  sigma <- sqrt(b_tilde/(a_tilde-1))
  print('Done!')
  mylist <- list(
    classic_coef = lm_classic$coef,
    classic_sigma = summary(lm_classic)$sigma,
    bayesian_coef = m_tilde,
    bayesian_sigma = sigma
  )
  return(mylist)
}
  
#' Bayesian Linear Regression with  noninformative prior.
#'
#' @param formula A string fromula.
#' @param data A data frame with your data including response.
#' @examples
#' test <- no_info(formula = 'dist ~ speed', data = cars)
#' test$classic_coef
#' test$classic_sigma
#' test$bayesian_coef
#' test$bayesian_sigma

no_info <- function(formula, data){
  lm_classic <- lm(formula, data)
  X <- model.matrix(lm_classic)
  y <- as.matrix(lm_classic$model[1])
  n <- dim(y)[1]
  p <- dim(X)[2]
  a <- -p
  M_tilde <- solve(t(X)%*%X)
  m_tilde <- M_tilde %*% t(X) %*% y
  a_tilde <- a + n/2
  b_tilde <- (1/2)*(t(y)%*%y - t(m_tilde)%*%t(X)%*%X%*%m_tilde)
  a_tilde <- as.numeric(a_tilde)
  b_tilde <- as.numeric(b_tilde)
  sigma <- sqrt(b_tilde/(a_tilde-1))
  print('Done!')
  mylist <- list(
    classic_coef = lm_classic$coef,
    classic_sigma = summary(lm_classic)$sigma,
    bayesian_coef = m_tilde,
    bayesian_sigma = sigma
  )
  return(mylist)
}
