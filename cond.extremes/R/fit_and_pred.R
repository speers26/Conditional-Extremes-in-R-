#
# Script containing functions for fitting H&T (2004) alongside any dependencies
#

# conflict resolutions ----------------------------------------------------
pfrechet <- evd::pfrechet
qfrechet <- evd::qfrechet
dfrechet <- evd::dfrechet
rfrechet <- evd::rfrechet


# functions ---------------------------------------------------------------

log_sample_lap <- function(sample_n,alpha){
  #' function which generates a bivariate sample with MEVD logistic dependence
  #' structure and laplace margins
  #'
  #' @param sample_n an integer, size of sample to be generated
  #' @param alpha floating point in [0,1], MEVD logistic dependence parameter
  #'
  #' @returns vector, sample with laplace margins and logistic dep structure
  #' @export
  sample <- rmvevd(sample_n,dep = alpha,d=2,model="log",mar=c(1,1,1)) ## with unit frechet margins
  l_sample <- qlaplace(pfrechet(sample)) ## with laplace margins
  return(l_sample)
}

ht_pred <- function(v, theta, Z, n_pred, Y=T, points=F, plot=T, col="darkred"){
  #' generates prediction points for a given H&T model
  #'
  #' @param v float, value of conditioned variable to predict above
  #' @param theta vector, pars from H&T fit to use in prediction
  #' @param Z vector, residuals from H&T fit to use in prediction
  #' @param n_pred int, number of predictions to make
  #' @param Y bool, true if we are modeling Y|X, false if X|Y
  #' @param plot bool, set to true to see visualisation
  #' @param col string, sets the color of plotted predictions
  #'
  #' @returns matrix of predicted X,Y points
  #' @export
  alpha <- theta[1] ; beta <- theta[2]

  rand_yi <- v + rexp(n_pred)
  rand_z <- sample(Z, n_pred, replace=T)

  bw = density(Z)$bw
  adjust = 0.1
  rand_z = rand_z + rnorm(length(rand_z), mean=0, sd= adjust*bw)

  prdctd_ymnsi <- rand_yi * alpha + rand_yi^beta * rand_z

  if(plot){
    if(1-Y){abline(h=v,col=col,lty="dashed")
      points(prdctd_ymnsi,rand_yi,col=col,pch=16,cex=0.5)}
    else{abline(v=v,col=col,lty="dashed")
      points(rand_yi,prdctd_ymnsi,col=col,pch=16,cex=0.5)}
  }

  if(points){
    return(matrix(data=c(rand_yi, prdctd_ymnsi), ncol=2))
  }
}

ht_fit = function(Yi, Ymnsi, qu, Y=T, keef=F, theta0, plot=T){
  #' fits the H&T (2004) model using Keef's (2013) constraints by optimising
  #' for alpha first, then beta
  #'
  #' @param Yi vector, data for conditioning variable
  #' @param Ymnsi vector, data for dependent variable
  #' @param q float, threshold quantile to fit H&T model above
  #' @param Y bool, true if we are modeling Y|X, false if X|Y
  #' @param theta0 vector of initial H&T parameter values to be used in optimisation
  #' @param keef bool, true if we are using Keef (2013) constraints
  #' @param plot bool, set to true to see visualisation
  #'
  #' @returns list, contains fitted parameters and residuals
  #' @export
  u <- quantile(Yi,qu)
  Ymnsi_rstr <- Ymnsi[Yi>=u] ; Yi_rstr <- Yi[Yi>=u]

  if(plot){ ## only need to run if just doing alpha first method
    if(1-Y){abline(h=u,lty="dashed",col="lightblue");
      points(Ymnsi_rstr,Yi_rstr,col="lightblue",pch=16,cex=0.5)}
    else{abline(v=u,lty="dashed",col="lightblue");
      points(Yi_rstr,Ymnsi_rstr,col="lightblue",pch=16,cex=0.5)}
  }

  if(1-keef){
    theta_beta0 <- optim(fn=b0_residualnLL,Ymnsi=Ymnsi_rstr,Yi=Yi_rstr,par=theta0, method="Nelder-Mead")$par
    theta_beta0 = c(theta_beta0[1], 0, theta_beta0[2], theta_beta0[3])
    #print(theta_beta0)
    init_theta = optim(fn=ht_residualnLL,Ymnsi=Ymnsi_rstr,Yi=Yi_rstr,par=theta_beta0, method="Nelder-Mead")$par
  }
  else{
    theta_beta0 <- optim(fn=b0_keef_residualnLL,Ymnsi=Ymnsi_rstr,Yi=Yi_rstr,par=theta0, method="Nelder-Mead")$par
    theta_beta0 = c(theta_beta0[1], 0, theta_beta0[2], theta_beta0[3])
    #print(theta_beta0)
    init_theta = optim(fn=ht_keef_residualnLL,Ymnsi=Ymnsi_rstr,Yi=Yi_rstr,par=theta_beta0, method="Nelder-Mead")$par
  }

  residuals <- (Ymnsi_rstr - init_theta[1] * Yi_rstr) / Yi_rstr^init_theta[2]
  return(list(theta=init_theta,res=residuals))
}


# internal functions ------------------------------------------------------
ht_residualnLL <- function(Ymnsi, Yi, theta){
  #' computes the negative residual log-likelihood for H&T model under G
  #' normality assumption
  #'@keywords internal
  #'
  #'@param Ymnsi vector for dependent variable
  #'@param Yi vector of conditioned variable
  #'@param theta vector of H&T parameters
  #'
  #'@returns float, negative log-likelihood
  alpha <- theta[1] ; beta <- theta[2] ; mu <- theta[3] ; sigma <- theta[4]
  if (beta < -1 | beta > 1 | alpha < -1 | alpha > 1 | sigma < 0){return(1e10)}
  else {
    LL <- sum(dnorm(Ymnsi,mean=alpha*Yi + (Yi^beta)*mu,sd=(Yi^beta)*sigma,log=T))
    return(-LL)
  }
}

ht_keef_residualnLL <- function(Ymnsi, Yi, theta){
  #' computes the negative residual log-likelihood for H&T model under G
  #' normality assumption with keef (2013) constraints
  #' @keywords internal
  #'
  #'@param Ymnsi vector for dependent variable
  #'@param Yi vector of conditioned variable
  #'@param theta vector of H&T parameters
  #'
  #'@returns float, negative log-likelihood
  alpha <- theta[1] ; beta <- theta[2] ; mu <- theta[3] ; sigma <- theta[4]
  cases_q0 <- both_cases_check(alpha, beta, Yi, Ymnsi, q = 1)
  cases_q1 <- both_cases_check(alpha, beta, Yi, Ymnsi, q = 0)
  if (cases_q0 && cases_q1){
    ### something causes NaNs to be produces here if v is large
    LL <- sum(dnorm(Ymnsi,mean=alpha*Yi + (Yi^beta)*mu,sd=(Yi^beta)*sigma,log=T))
    ###
    return(-LL)
  }
  else{return(1e10)}
}

b0_residualnLL = function(Ymnsi, Yi, theta){ ## here theta is 3d
  #' computes the negative residual log-likelihood for H&T model under G
  #' normality assumption, with beta=0 fixed
  #' @keywords internal
  #'
  #'@param Ymnsi vector for dependent variable
  #'@param Yi vector of conditioned variable
  #'@param theta vector of H&T parameters
  #'
  #'@returns float, negative log-likelihood
  alpha <- theta[1] ; beta <- 0 ; mu <- theta[2] ; sigma <- theta[3]
  if (alpha < -1 | alpha > 1 | sigma <= 0){return(1e10)}
  else {
    LL <- sum(dnorm(Ymnsi,mean=alpha*Yi + (Yi^beta)*mu,sd=(Yi^beta)*sigma,log=T))
    return(-LL)
  }
}

b0_keef_residualnLL = function(Ymnsi, Yi, theta){
  #' computes the negative residual log-likelihood for H&T model under G
  #' normality assumption, with beta=0 fixed and keef (2013) constraints
  #' @keywords internal
  #'
  #'@param Ymnsi vector for dependent variable
  #'@param Yi vector of conditioned variable
  #'@param theta vector of H&T parameters
  #'
  #'@returns float, negative log-likelihood
  alpha <- theta[1] ; beta <- 0 ; mu <- theta[2] ; sigma <- theta[3]
  cases_q0 <- both_cases_check(alpha, beta, Yi, Ymnsi, q = 1)
  cases_q1 <- both_cases_check(alpha, beta, Yi, Ymnsi, q = 0)
  if (cases_q0 && cases_q1){
    ### something causes NaNs to be produces here if v is large
    LL <- sum(dnorm(Ymnsi,mean=alpha*Yi + (Yi^beta)*mu,sd=(Yi^beta)*sigma,log=T))
    ### -------------------------------------------------------
    return(-LL)
  }
  else{return(1e10)}
}

both_cases_check <- function(alpha, beta, yi, ymnsi, q){
  #' computes required quantiles and feeds them to case one  adn two functions
  #' given by Keef (2013)
  #' @keywords internal
  #'
  #' @param alpha float, H&T alpha parameter
  #' @param beta float, H&T beta parameter
  #' @param yi vector of conditioning variable
  #' @param ymnsi vector of dependent variable
  #' @param q ???? check Keef (2013) for this
  #'
  #' @returns boolean, true if case on is satisfied
  ## get quantiles
  # old constraints
  if( alpha < -1 | alpha > 1 | beta < -1 | beta > 1){return(0)}

  ## get quantiles
  z_pos <- ymnsi - yi
  z_neg <- ymnsi + yi
  z_ind <- (ymnsi - alpha * yi)/yi^beta

  q_pos <- quantile(z_pos, q, names=F)
  q_neg <- quantile(z_neg, q, names=F)
  q_ind <- quantile(z_ind, q, names=F)

  quantiles <- list(pos = q_pos, neg = q_neg, ind = q_ind)

  theta <- list(alpha = alpha, beta = beta)

  v <- max(yi) + 10 ## increasing this seems to give NaN error sometimes in dnorm

  ## check new constraints
  return(case_one(theta, quantiles, v)*case_two(theta, quantiles, v))
}

case_one <- function(theta, quantiles, v){
  #' returns truth value of case one for constraints given by keef (2013)
  #' @keywords internal
  #'
  #' @param theta H&T parameter
  #' @param quantiles ?? check Keef(2013) for this
  #' @param v ?? check keef(2013) for this
  #'
  #' @returns boolean, true if case 1 satisfied
  alpha <- theta$alpha ; beta <- theta$beta

  q_pos <- quantiles$pos
  q_ind <- quantiles$ind

  term1 <- 1
  term2 <- 1 - beta * q_ind * v^(beta - 1)
  term3 <- 1 - v^(beta - 1) * q_ind + v^-1 * q_pos

  cond1 <- alpha <= min(term1, term2, term3)

  term4 <- 1 - beta * q_ind^(beta - 1)
  term5 <- (1 - beta^-1) * (beta * q_ind)^(1/(1-beta)) * (1 - alpha)^(-beta/(1-beta)) + q_pos

  cond2 <- (term4 < alpha && alpha <= 1) * (term5 > 0)

  ## if cannt evaluate condition take as false
  if(is.na(cond1)){cond1<-0}
  if(is.na(cond2)){cond2<-0}

  return(cond1 + cond2)

}

case_two <- function(theta, quantiles, v){
  #' returns truth value of case two for constraints given by keef (2013)
  #' @keywords internal
  #'
  #' @param theta H&T parameter
  #' @param quantiles ?? check Keef(2013) for this
  #' @param v ?? check keef(2013) for this
  #'
  #' @returns boolean, true if case 2 satisfied
  alpha <- theta$alpha ; beta <- theta$beta

  q_neg <- quantiles$neg
  q_ind <- quantiles$ind

  term1 <- 1
  term2 <- 1 + beta * v^(beta - 1) * q_ind
  term3 <- 1 + v^(beta - 1) * q_ind - v^-1 * q_neg

  cond1 <- -alpha <= min(term1, term2, term3)

  term4 <- 1 + beta * v^(beta - 1) * q_ind
  term5 <- (1 - beta^-1) * (-beta * q_ind)^(1/(1-beta)) * (1 + alpha)^(-beta/(1-beta)) - q_neg

  cond2 <- (term4 < -alpha && -alpha <= 1) * (term5 > 0)

  ## if cannt evaluate condition take as false
  if(is.na(cond1)){cond1<-0}
  if(is.na(cond2)){cond2<-0}

  return(cond1 + cond2)

}
