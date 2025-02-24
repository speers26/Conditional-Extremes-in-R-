#
# Function for using H&T (2004)/ Keef(2013) to estimate extreme joint exceedance probabilities
#

prob_est = function(data, v, q.cond, npred=1e5, alt=FALSE){
  #' Probability estimation using the conditional extremes model
  #' 
  #' @param data a matrix with two columns, the first being the x values and the second being the y values, on laplace margins
  #' @param v the joint exceedance value on laplace margins
  #' @param q.cond vector of conditioning quantiles
  #' @param npred number of predictions to make in empirical probability estimation
  #' @param alt logical, if TRUE, the function will use the alternative method for estimating conditional extremes parameters
  #' 
  #' @return the estimated probability of joint exceedance 
  #' @export
  
  exp_sample = rexp(npred)
  
  qu_x = q.cond[1]
  qu_y = q.cond[2]
  
  if (alt){
    
    x2condx1 = ht.fitAlt(data[,1], data[,2], qu_x)
    x1condx2 = ht.fitAlt(data[,2], data[,1], qu_y)
    
  } else {
    
    x2condx1 = ht.fit(data[,1], data[,2], qu_x, keef=TRUE, theta0=c(0,0,1), plot=FALSE)
    x1condx2 = ht.fit(data[,2], data[,1], qu_y, keef=TRUE, theta0=c(0,0,1), plot=FALSE)
    
  }
 
  z1 = sample(x2condx1$res, npred, replace=TRUE) 
  z2 = sample(x1condx2$res, npred, replace=TRUE) 
  
  x2condx1_xpred = exp_sample + v
  x2condx1_ypred = x2condx1$theta[1] * x2condx1_xpred + z1 * x2condx1_xpred^x2condx1$theta[2]
  
  x2condx1_prob = sum((x2condx1_ypred > v) & (x2condx1_xpred > x2condx1_ypred)) / length(z1)
  
  x1condx2_ypred = exp_sample + v
  x1condx2_xpred = x1condx2$theta[1] * x1condx2_ypred + z2 * x1condx2_ypred^x1condx2$theta[2]
  x1condx2_prob = sum((x1condx2_xpred > v) & (x1condx2_ypred > x1condx2_xpred)) / length(z1)
  
  joint_prob = (x2condx1_prob + x1condx2_prob) * 0.5 * exp(-v) 
  
  return(joint_prob)
  
}
