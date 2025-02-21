prob_est = function(data, v, q.cond, npred=1e5){
  #' Probability estimation using the conditional extremes model
  #' 
  #' @param data a matrix with two columns, the first being the x values and the second being the y values, on laplace margins
  #' @param v the joint exceedance value on laplace margins
  #' @param q.cond vector of conditioning quantiles
  #' @param npred number of predictions to make in empirical probability estimation
  #' 
  #' @return the estimated probability of joint exceedance 
  #' @export
  
  # fit conditional extremes models
  exp_sample = rexp(npred)
  
  qu_x = q.cond[1]
  qu_y = q.cond[2]
  
  x2condx1 = ht.fit(data[,1], data[,2], qu_x, keef=TRUE, theta0=c(0,0,1), plot=FALSE)
  x1condx2 = ht.fit(data[,2], data[,1], qu_y, keef=TRUE, theta0=c(0,0,1), plot=FALSE)
  z1 = sample(x2condx1$res, npred, replace=TRUE) # + rnorm(npred, sd=0.1)
  z2 = sample(x1condx2$res, npred, replace=TRUE) # + rnorm(npred, sd=0.1)
  
  x2condx1_xpred = exp_sample + v
  x2condx1_ypred = x2condx1$theta[1] * x2condx1_xpred + z1 * x2condx1_xpred^x2condx1$theta[2]
  
  x2condx1_prob = sum((x2condx1_ypred > v) & (x2condx1_xpred > x2condx1_ypred)) / length(z1)
  
  x1condx2_ypred = exp_sample + v
  x1condx2_xpred = x1condx2$theta[1] * x1condx2_ypred + z2 * x1condx2_ypred^x1condx2$theta[2]
  x1condx2_prob = sum((x1condx2_xpred > v) & (x1condx2_ypred > x1condx2_xpred)) / length(z1)
  
  joint_prob = (x2condx1_prob + x1condx2_prob) * 0.5 * exp(-v) 
  
  return(joint_prob)
  
}
