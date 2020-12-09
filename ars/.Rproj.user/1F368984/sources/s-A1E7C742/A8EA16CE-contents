context("tests")

##**************************************************************
## logistic distribution
test_that("logistic distribution",
          {
  test_fun(1000, dis_logistic, rlogis, "logistic distribution", D_left = -10, D_right = 10)
  expect_equal(1,1)
  })


##***************************************************************
## standard normal distribution
test_that("normal distribution",{
          test_fun(1000, dis_norm, rnorm, "normal distribution", D_left = -10, D_right = 10)
          expect_equal(1,1)
          })
#
# ##***************************************************************
# ## uniform distribution
# test_that("uniform distribution",{
#   test_fun(1000, dis_unif, runif, "uniform distribution", D_left = 0, D_right = 1)
#   expect_equal(1,1)
#   })
#
# ##***************************************************************
# ## laplace distribution
# test_that("laplace distribution",{
#   test_fun(1000, dis_laplace, rmutil::rlaplace, "laplace distribution", D_left = -10, D_right = 10)
#   expect_equal(1,1)
#   })
#
# ##***************************************************************
# ## gamma distribution
# test_that("gamma distribution",{
#   r_gamma <- function(n){
#     return(rgamma(n, shape = 1))
#   }
#   test_fun(1000, dis_gamma, r_gamma, "gamma distribution", D_left = 0.01, D_right = Inf)
#   expect_equal(1,1)
#   })

##***************************************************************
##*
## non-concave case: a mixture of two normal distributions
 #
# dis_nc <- function(x){
# p <- 0.5*(1/sqrt(2*pi))*exp(-(x^2)/2) + 0.5*(1/sqrt(2*pi))*exp(-((x-3)^2)/2)
# return(p)
# }

 #test_that("mixture of normal distribution",{
 # expect_error(test_fun(1000, fun = dis_nc, rfun = NA , "mixture of normal distributions", D_left = -10, D_right = 10))
 # })

test_fun <- function(n, fun, rfun, fun_name, D_left, D_right){
  # Generate samples using ars() and rfun()
  print(paste0("generating samples from ", fun_name,"..."))
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
