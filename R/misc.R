##
## script containing extra functions which are useful when fitting H&T(2004)
##

pspliced = function(x0, x, u, gpd_par){
  #' Distribution function for GPD spliced with empirical
  #'
  #' spliced distribution which uses empirical cdf for below threshold and gpd
  #' for above
  #'
  #' @param x0 float to evaluate the distribution at
  #' @param x vector of sample used to fit the distribution
  #' @param u float, exceedance threshold
  #' @param gpd_par vector of gpd scale and shape parameters
  #'
  #' @returns float between [0,1], estimated cdf at x0
  #' @export
  if(x0 <= u){
    return(sum(x<=x0)/length(x))
  }
  else{
    p_exceed = 1 - pspliced(u, x, u, gpd_par)
    excess = x0 - u
    gpd_surv = evd::pgpd(excess, scale=gpd_par[1], shape=gpd_par[2], lower.tail=F)
    return(1 - (gpd_surv * p_exceed))
  }
}


qspliced = function(p, x, q, gpd_par){
  #' Quantile function for GPD spliced with empirical
  #'
  #' spliced inverse distribution which uses empirical cdf for below threshold
  #' and gpd for above
  #'
  #' @param p float between [0,1] for which to find the corresponding quantile
  #' @param x vector of sample used to fit the distribution
  #' @param q float between [0,1], exceedance quantile
  #' @param gpd_par vector of gpd scale and shape parameters
  #'
  #' @returns float, estimated pth quantile for data x with gpd tail
  #' @export
  if(p <= q){
    return(quantile(x, p, names=F))
  }
  else{
    u = quantile(x, q, names=F)
    gpd_prob = 1 - (1-p)/(1-q)
    gpd_quant = evd::qgpd(gpd_prob, scale=gpd_par[1], shape=gpd_par[2])
    return(u + gpd_quant)
  }
}
