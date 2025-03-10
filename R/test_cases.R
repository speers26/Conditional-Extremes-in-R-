#
# Sampling from representative bivariate distributions with Laplace margins.
#


# conflict resolutions ----------------------------------------------------
pfrechet <- evd::pfrechet
qfrechet <- evd::qfrechet
dfrechet <- evd::dfrechet
rfrechet <- evd::rfrechet

# functions for sample generation  ---------------------------------------------------------------
bivariate_gaussian_laplace = function(sample_n, rho){
  #' Bivariate Gaussian with Laplace margins
  #' 
  #' function which generates a bivariate Gaussian sample with laplace margins
  #' 
  #' @param sample_n an integer, size of sample to be generated
  #' @param rho float, correlation between the two variables
  #' 
  #' @returns matrix, sample with laplace margins and gaussian dependence structure
  #' @export
  
  # simulate bivariate Gaussian data
  mu = c(0, 0)
  Sigma = matrix(c(1, rho, rho, 1), nrow = 2)
  X = rmvnorm(n = sample_n, mean = mu, sigma = Sigma)
  
  # change to laplace margins
  X = qlaplace(pnorm(X))
  
  return(X)
  
}

neglog_sample_lap = function (sample_n, alpha){
  #' MEVT sample with laplace margins generation
  #' 
  #' function which generates a bivariate sample with MEVD inverted logistic dependence
  #' structure and laplace margins
  #' 
  #' @param sample_n an integer, size of sample to be generated
  #' @param alpha floating point in [0,1], MEVD inverted logistic dependence parameter
  #' 
  #' @returns matrix, sample with laplace margins and inverted logistic dep structure
  #' @export

  sample <- rmvevd(sample_n, dep = alpha, d = 2, model = "log", 
                   mar = c(1, 1, 1))
  l_sample <- -qlaplace(evd::pfrechet(sample))
  return(l_sample)
}


log_sample_lap <- function(sample_n, alpha){
  #' MEVT sample with laplace margins generation
  #'
  #' function which generates a bivariate sample with MEVD logistic dependence
  #' structure and laplace margins
  #'
  #' @param sample_n an integer, size of sample to be generated
  #' @param alpha floating point in [0,1], MEVD logistic dependence parameter
  #'
  #' @returns matrix, sample with laplace margins and logistic dep structure
  #' @export
  
  sample <- rmvevd(sample_n,dep = alpha,d=2,model="log",mar=c(1,1,1)) ## with unit frechet margins
  l_sample <- qlaplace(pfrechet(sample)) ## with laplace margins
  return(l_sample)
}

# functions for true joint exceedance probabilities -----------------------
logistic_prob = function(v, alpha){
  #' true probability of joint exceedance for a logistic dependence structure
  #' 
  #' @param alpha the dependence parameter
  #' @param v the joint exceedance value on laplace margins
  #' 
  #' @return the true probability of joint exceedance
  #' @export
  
  v.frechet = evd::qfrechet(plaplace(v))
  
  p = 1 - 2 * evd::pfrechet(v.frechet) + pmvevd(q=c(v.frechet, v.frechet), dep=alpha,
                                                model="log", d=2, mar=c(1,1,1), lower.tail=TRUE)
  return(p)
    
}

gaussian_prob = function(v, rho){
  #' true probability of joint exceedance for a gaussian dependence structure
  #' 
  #' @param rho the dependence parameter
  #' @param v the joint exceedance value on laplace margins
  #' 
  #' @return the true probability of joint exceedance
  #' @export
  
  v.gaussian = v %>% plaplace %>% qnorm
  
  p = 1 - 2 * pnorm(v.gaussian) + pmvnorm(upper=c(v.gaussian, v.gaussian), mean=c(0,0), corr=matrix(data=c(1, rho, rho, 1), nrow=2))[[1]]
  
  return(p)
  
} 

neglog_prob = function(v, alpha){
  #' true probability of joint exceedance for an inverted logistic dependence structure
  #' 
  #' @param alpha the dependence parameter
  #' @param v the joint exceedance value on laplace margins
  #'
  #' @return the true probability of joint exceedance
  #' @export
  
  v.frechet = evd::qfrechet(plaplace(-v))
  
  p = pmvevd(q=c(v.frechet, v.frechet), dep=alpha, model="log", d=2, mar=c(1,1,1), lower.tail=TRUE)
  
  return(p)

}
  