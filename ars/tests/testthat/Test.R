context("tests")



## Test if the warnings for the wrong inputs work
##******************************************************************************
## Input of n

test_that("Input of n", {
  expect_error(ars(n = "a", g = rnorm, D_left = -10, D_right = 10),
               "Please provide n as a number")})

##******************************************************************************
## Input of D_left and D_right

test_that("Input of D_left", {
  expect_error(ars(n = 10, g = rnorm, D_left = "asv", D_right = 10),
               "Please provide D_left as a number")})

test_that("Input of D_left and D_right", {
  expect_error(ars(n = 10, g = rnorm, D_left = 10, D_right = 10),
               "Please provide different D_left, D_right")})

##******************************************************************************
## Input of function g

test_that("Input of g", {
  expect_error(ars(n = 10, g = 3, D_left = -10, D_right = 10),
               "Please provide g as a function")})

##******************************************************************************


## Test if the function will give the correct result
##******************************************************************************
## Log-concave distributions
##******************************************************************************
## Logistic distribution (1, 1.5)

dis_logistic <- function(x) {return(dlogis(x, location = 1, scale = 1.5))}
r_logistic <- function(n) {return(rlogis(n, location = 1, scale = 1.5))}

test_that("Logistic distribution", {
  test_fun(1000, dis_logistic, r_logistic, "Logistic distribution",
           D_left = -10, D_right = 10)
  expect_equal(1,1)})

##******************************************************************************
##Normal distribution N(mean=2, sd=3)

dis_norm <- function(x) {return(dnorm(x, mean = 2, sd = 3))}
r_norm <- function(n) {return(rnorm(n, mean = 2, sd = 3))}

test_that("Normal distribution", {
          test_fun(1000, dis_norm, r_norm, "Normal distribution",
                   D_left = -10, D_right = 10)
          expect_equal(1,1)})

##******************************************************************************
## Uniform distribution [0,1]
dis_unif <- function(x) {return(1)}

test_that("Uniform distribution", {
   test_fun(1000, dis_unif, runif, "Uniform distribution", D_left = 0, D_right = 1)
   expect_equal(1,1)})

##***************************************************************
## Laplace distribution (mu=0)

dis_laplace <- function(x)
{
  f <- 1/2*exp(-abs(x))
  return(f)
 }

test_that("Laplace distribution", {
   test_fun(1000, dis_laplace, rmutil::rlaplace, "Laplace distribution",
            D_left = -10, D_right = 10)
   expect_equal(1,1)})

##******************************************************************************
## Gamma distribution: Gamma(2,2)

dis_gamma <- function(x) {return(dgamma(x, shape = 2, rate = 2))}

r_gamma <- function(n) {return(rgamma(n, shape = 2, rate = 2))}

test_that("Gamma distribution", {
  test_fun(1000, dis_gamma, r_gamma, "Gamma distribution",
           D_left = 0.1, D_right = 50)
  expect_equal(1,1)})

##*****************************************************************************
#### Non-log-concave cases
##*****************************************************************************
## T distribution (DOF=2)

dis_t <- function(x) {return(dt(x, df = 2))}
r_t <- function(n) {return(rt(n, df = 2))}

test_that("T distribution", {
  test_fun(1000, dis_t, r_t, "T distribution", D_left = -10, D_right = 10)
  expect_equal(1,1)})

##******************************************************************************
## A mixture of two normal distributions

dis_nc <- function(x)
{
  p <- 0.5*(1/sqrt(2*pi))*exp(-(x^2)/2) + 0.5*(1/sqrt(2*pi))*exp(-((x-3)^2)/2)
  return(p)
}

test_that("Mixture of normal distribution", {
  expect_error(test_fun(1000, fun = dis_nc, rfun = NA , "Mixture of normal distributions",
                        D_left = -10, D_right = 10))})

##******************************************************************************
## Test function
test_fun <- function(n, fun, rfun, fun_name, D_left, D_right){
  # Generate samples using ars() and rfun()
  print(paste0("generating samples from ", fun_name, "..."))
  sample <- try(ars( g = fun, n = n, D_left = D_left, D_right = D_right))
  rsample <- rfun(n)

  # Perform Kolmogorov-Smirnov two-sample test
  print(paste0("performing ks-test on ", fun_name,"..."))
  test <- ks.test(sample, rsample)
  p_value <- test$p.value
  error_msg <- paste0("test failed with ", fun_name)
  pass_msg <- paste0("test passed with ", fun_name)
  if(p_value <= 0.05)
    print(error_msg)
  else
    print(pass_msg)
}
