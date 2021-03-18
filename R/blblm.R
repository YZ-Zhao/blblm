#' @title Bootstrap for linear regression model
#' @description Bag of little bootstraps for linear regression model.

#' @import purrr
#' @import stats
#' @import furrr
#' @import future
#' @import parallel
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
#'
#'

"_PACKAGE"



## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' "blblm()" is a function to model linear regression with Little Bag of Bootstraps
#' @param formula A symbolic description of the model
#' @param data dataframe
#' @param m integer,number of subsample you want
#' @param B integer, bootstrap iteration times
#' @param nthreads the number of cores you want to use when processing
#' @examples
#' blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
#' @return object
#' @export
blblm <- function(formula, data, m = 10, B = 5000, nthreads = 1) {
  data_list <- split_data(data, m)
  if(nthreads > 1){
    suppressWarnings(plan(multiprocess, workers = nthreads))
    options(future.rng.onMisuse = "ignore")
    estimates <- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  }
  else{
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  }
  if(!inherits(plan(), "sequential")) plan(sequential)
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#' "split_data()" split data into m parts of approximated equal sizes
#' The sampling is with replacement.
#' @param data dataframe
#' @param m integer, number of dataframes you want to split to
#' @return a list of dataframes
#'
#' @export
#'
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}



#' Linear regression modeling based on one subsample
#' @param formula A symbolic description of the model
#' @param data dataframe
#' @param n sample size
#' @param B integer, bootstrap iteration times
lm_each_subsample <- function(formula, data, n, B) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  m <- model.frame(formula, data)
  X <- model.matrix(formula, m)
  y <- model.response(m)
  replicate(B, lm1(X, y, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
#'
#' @param X A symbolic description of the model
#' @param y dataframe
#' @param n freqency of each observations
#'
#' @return list, coefficients and sigma
#'
#' @export
lm1 <- function(X, y, n) {
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  fit <- lm.wfit(X, y, freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
#' @param fit object of blblm
#' @return coefficient of the model
#' @export
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#' @param fit fitted regression
#'
#' @return estimated standard deviation of the fitted model
#'
blbsigma <- function(fit) {
  p <- fit$rank
  e <- fit$residuals
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' Print the symbolic model
#' @param x object
#'
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}

#' Get the sd of blblm
#' @param object an object of "blblm"
#' @param confidence logic, need confidence interval or not
#' @param level the confidence level
#' @param nthreads the number of cores you want to use when processing
#' @return When confidence = FALSE, return estiamted sd, else return sd and its upper bound and lower bound
#'
#' @export

#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, nthreads = 1) {
    est <- object$estimates
    sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
    if(nthreads > 1){
      suppressWarnings(plan(multiprocess, workers = nthreads))
      options(future.rng.onMisuse = "ignore")
      if (confidence) {
        alpha <- 1 - 0.95
        limits <- est %>%
          future_map_mean(~ quantile(future_map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
          set_names(NULL)
        return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
      } else {
        return(sigma)
      }
    }
    else{
      if (confidence) {
        alpha <- 1 - 0.95
        limits <- est %>%
          map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
          set_names(NULL)
        return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
      } else {
        return(sigma)
    }
    }
}



#' Get the coefficients
#' @param object blblm object
#'
#' @return coefficients
#'
#' @export
#' @method coef blblm

coef.blblm <- function(object) {
    est <- object$estimates
    map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
  }


#' obtain confidence interval of coefficients
#'
#' @param object object of "blblm"
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level confidence level, the default value is 0.95
#' @param nthreads the number of cores you want to use when processing
#' @return the confidence intervals of coefficients
#'
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, nthreads = 1) {
    if (is.null(parm)) {
      parm <- attr(terms(object$formula), "term.labels")
    }
    alpha <- 1 - level
    est <- object$estimates
    if(nthreads > 1){
      suppressWarnings(plan(multiprocess, workers = nthreads))
      options(future.rng.onMisuse = "ignore")
      out <- future_map_rbind(parm, function(p) {
        future_map_mean(est, ~ future_map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
      })
    }
    else{
      out <- map_rbind(parm, function(p) {
        map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
      })
    }
    if (is.vector(out)) {
      out <- as.matrix(t(out))
    }
    dimnames(out)[[1]] <- parm
    out
}

#' prediciton
#'
#' @param object an object of "blblm"
#' @param new_data new data frame
#' @param confidence need confidence interval or not
#' @param level confidence level, default is 0.95
#' @param nthreads the number of cores you want to use when processing#' @param nthreads the number of cores you want to use when processing

#' @return prediction result
#'
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, nthreads = 1) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if(nthreads > 1){
    suppressWarnings(plan(multiprocess, workers = nthreads))
    options(future.rng.onMisuse = "ignore")
    if (confidence) {
      future_map_mean(est, ~ future_map_cbind(., ~ X %*% .$coef) %>%
        apply(1, mean_lwr_upr, level = level) %>%
        t())
    } else {
      future_map_mean(est, ~ future_map_cbind(., ~ X %*% .$coef) %>% rowMeans())
    }
  }
  else{
    if (confidence) {
     map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
                        apply(1, mean_lwr_upr, level = level) %>%
                        t())
    } else {
      map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
    }
  }
}



mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

future_map_mean <- function(.x, .f, ...) {
  (future_map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}



map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

future_map_cbind <- function(.x, .f, ...) {
  future_map(.x, .f, ...) %>% reduce(cbind)
}


map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}

future_map_rbind <- function(.x, .f, ...) {
  future_map(.x, .f, ...) %>% reduce(rbind)
}


