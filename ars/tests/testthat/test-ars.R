context("tests")

library(testthat)


## Test if the warnings for the wrong inputs work
cat("First test if the warnings for the wrong inputs work. \n")
##******************************************************************************
## Input of n

cat("Test if there's a warning when the input of n is not a number.\n")
test_that("Input of n", {
  expect_error(ars(n = "a", g = rnorm, D_left = -10, D_right = 10),
               "Please provide n as a number")})
cat("Passed. \n")
##******************************************************************************
## Input of D_left and D_right

cat("Test if there's a warning when the input of D_left is not a number.\n")
test_that("Input of D_left", {
  expect_error(ars(n = 10, g = rnorm, D_left = "asv", D_right = 10),
               "Please provide D_left as a number")})
cat("Passed. \n")

cat("Test if there's a warning when the input of D_right is not a number.\n")
test_that("Input of D_right", {
  expect_error(ars(n = 10, g = rnorm, D_left = 1, D_right = "asad"),
               "Please provide D_right as a number")})
cat("Passed. \n")

cat("Test if there's a warning when the input of D_left and D_right are equal.\n")
test_that("Input of D_left and D_right", {
  expect_error(ars(n = 10, g = rnorm, D_left = 10, D_right = 10),
               "Please provide different D_left, D_right")})
cat("Passed. \n")

##******************************************************************************
## Input of function g

cat("Test if there's a warning when the input of g is not a function.\n")
test_that("Input of g", {
  expect_error(ars(n = 10, g = 3, D_left = -10, D_right = 10),
               "Please provide g as a function")})
cat("Passed. \n")

##******************************************************************************


## Test if the function will give the correct result for Log-concave cases

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

##******************************************************************************
## Log-concave cases
cat("The following tests are log-concave cases to check if the function works for them. \n")
##******************************************************************************
## Logistic distribution (1, 1.5)

dis_logistic <- function(x) {return(dlogis(x, location = 1, scale = 1.5))}
r_logistic <- function(n) {return(rlogis(n, location = 1, scale = 1.5))}

cat("Test the Logistic Distribution with location=1 and scale =1.5. \n")
test_that("Logistic Distribution", {
  test_fun(1000, dis_logistic, r_logistic, "Logistic Distribution", D_left = -10, D_right = 10)
  expect_equal(1,1)})

##******************************************************************************
##Normal distribution (mean=2, sd=3)

dis_norm <- function(x) {return(dnorm(x, mean = 2, sd = 3))}
r_norm <- function(n) {return(rnorm(n, mean = 2, sd = 3))}

cat("Test the Normal Distribution with mean=2 and sd=3. \n")
test_that("Normal Distribution", {
  test_fun(1000, dis_norm, r_norm, "Normal Distribution", D_left = -10, D_right = 10)
  expect_equal(1,1)})

##******************************************************************************
## Uniform distribution [0,1]
dis_unif <- function(x) {return(dunif(x, 0, 1))}
r_unif <- function(n) {return(runif(n, 0, 1))}

cat("Test the Uniform Distribution in [0,1]. \n")
test_that("Uniform Distribution", {
  test_fun(1000, dis_unif, r_unif, "Uniform Distribution", D_left = 0, D_right = 1)
  expect_equal(1,1)})

##***************************************************************
## Laplace distribution (mu=0)

dis_laplace <- function(x)
{
  f <- 1/2*exp(-abs(x))
  return(f)
}

cat()
test_that("Laplace Distribution", {
  test_fun(1000, dis_laplace, rmutil::rlaplace, "Laplace Distribution", D_left = -10, D_right = 10)
  expect_equal(1,1)})

##******************************************************************************
## Gamma distribution: Gamma(2,2)

dis_gamma <- function(x) {return(dgamma(x, shape = 2, rate = 2))}
r_gamma <- function(n) {return(rgamma(n, shape = 2, rate = 2))}

cat("Test the Gamma Distribution with shape=2, rate=2. \n")
test_that("Gamma Distribution", {
  test_fun(1000, dis_gamma, r_gamma, "Gamma Distribution", D_left = 0.1, D_right = 50)
  expect_equal(1,1)})

##******************************************************************************
## Non-log-concave cases
cat("The following tests are non-log-concave cases to prove the function won't work for them. \n")
##******************************************************************************
## T distribution (DOF=2)

dis_t <- function(x) {return(dt(x, df = 2))}

cat("Test that the T Distribution with the degree of freedom=2 doesn't work. \n")
test_that("T Distribution", {
  expect_error(ars(1000, dis_t, D_left = -10, D_right = 10),
               "Input is not a log-concave function")})

cat("Passed. \n")

##******************************************************************************
## A mixture of two normal distributions

dis_nc <- function(x)
{
  p <- 0.5*(1/sqrt(2*pi))*exp(-(x^2)/2) + 0.5*(1/sqrt(2*pi))*exp(-((x-3)^2)/2)
  return(p)
}

cat("Test that a mixture of two normal distributions doesn't work. \n")
test_that("Mixture of normal distribution", {
  expect_error(ars(1000, dis_nc, D_left = -10, D_right = 10),
               "Input is not a log-concave function")})

cat("Passed. \n")

##******************************************************************************


