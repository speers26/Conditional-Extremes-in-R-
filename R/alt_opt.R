# alternate method --------------------------------------------------------

ht.fitAlt = function(Yi, Ymnsi, qu){
  #' alternative optimisation for H&T (2004) model over the constrained parameter space of Keef (2013)
  #' 
  #' optimisating in this way is faster, useful for cases where many repeated fits are needed (e.g., parametric bootstrapping)
  #' 
  #' @param Yi the excesses of the conditioning variable
  #' @param Ymnsi the excesses of the conditioned variable
  #' @param qu the quantile of the conditioning variable
  #' 
  #' @return list of the estimated parameters and residuals
  #' @export
  
  u = quantile(Yi, qu, names=FALSE)
  excess_x = Yi[Yi>=u] 
  excess_y = Ymnsi[Yi>=u]
  
  unCnstr_fit = optim(fn=ht_residualnLL, par=c(0.5,0.5,0,1), Yi=excess_x, Ymnsi=excess_y)
  unCnstr_theta = unCnstr_fit$par
  alpha_opt = unCnstr_theta[1]
  
  cases_q0 <- both_cases_check(unCnstr_theta[1], unCnstr_theta[2], excess_x, excess_y, q = 1)
  cases_q1 <- both_cases_check(unCnstr_theta[2], unCnstr_theta[2], excess_x, excess_y, q = 0)
  
  if (!cases_q0 | !cases_q1){
    
    message("Outside Keef constraints, finding closest point in parameter space")
    
    beta_max = unCnstr_theta[2]
    
    beta1 = optimise(equality1, c(0, beta_max), Yi=excess_x, Ymnsi=excess_y, alpha=unCnstr_theta[1], v=max(Yi) + 10)$minimum
    beta2 = optimise(equality2, c(0, beta_max), Yi=excess_x, Ymnsi=excess_y, alpha=unCnstr_theta[1], v=max(Yi) + 10)$minimum
    
    beta_opt = min(beta1, beta2)
    
  } else {
    
    beta_opt = unCnstr_theta[2]
    
  }
  
  res = (excess_y - alpha_opt*excess_x) / (excess_x^beta_opt)
  
  return(list(theta = c(alpha_opt, beta_opt), res=res))
  
}