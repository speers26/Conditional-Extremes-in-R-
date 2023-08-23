#
# Script containing functions for selecting threshold for H&T (2004)
#

# conflict resolutions ----------------------------------------------------
pfrechet <- evd::pfrechet
qfrechet <- evd::qfrechet
dfrechet <- evd::dfrechet
rfrechet <- evd::rfrechet

# functions ---------------------------------------------------------------

energy.thrsh.test = function(data, qt.range, n.boot, n.boot.par, n.rep){
  #' Energy independence test for threshold selection
  #'
  #' conducts test of independence based on that of Wan (2019) to validate H&T
  #' exceedance quantiles
  #' @param data matrix of 2 columns containing to fit HT to
  #' @param qt.range range of exceedance threshold quantiles to evaluate
  #' @param n.boot size of bootstraps to take at each rep
  #' @param n.boot.par size of boostrapped parameter samples to calculate at each rep
  #' @param n.rep number of p value calculation reps per quantiles
  #' @returns returns nothing
  #' @export
  #'
  # set up empty arrays
  pvs = matrix(data=0, ncol=length(qt.range), nrow=n.rep)
  mean.pvs = c()

  # iterate over quantiles and reps
  for (i in 1:length(qt.range)){
    print(i)
    for (j in 1:n.rep){
      # get bootstrapped X and Y
      rand.ind = sample(1:length(data[,1]), n.boot)
      rand.X = data[,1][rand.ind]
      rand.Y = data[,2][rand.ind]

      # get bootstrapped ht parameters
      par.boot = ht.par.boot(data=matrix(data=c(rand.X, rand.Y), ncol=2),
                             q=qt.range[i], n.boot=n.boot.par)

      # get excesses
      u = quantile(rand.X, qt.range[i], names=F)
      exs = rand.X[rand.X>u] - u

      # get 'average' residuals
      alpha.mean = mean(sapply(par.boot$theta, '[[', 1))
      beta.mean = mean(sapply(par.boot$theta, '[[', 2))
      res = (rand.Y[rand.X>u] - alpha.mean*exs) / exs^beta.mean

      # transform margins
      q.exs = cdf(exs);q.exs[q.exs==1] = 1-1e-10;q.exs[q.exs==0] = 1e-10
      q.res = cdf(res);q.res[q.res==1] = 1-1e-10;q.res[q.res==0] = 1e-10
      exs = qexp(q.exs)
      res = qexp(q.res)

      # do independence energy-test
      test = indep.test(exs, res, R=1000)
      pvs[j, i] = test$p.value
    }
  }
  # find mean p values
  mean.pvs = apply(pvs, 2, mean)

  # plot mean p values
  plot(rev(1-qt.range), rev(mean.pvs), pch=2, col="red",
       ylim=c(min(pvs), max(pvs)), xlab="qt", ylab="p-values")

  # plot all p values
  for (q in rev(1:length(qt.range))){
    points(rep(1-qt.range[q], n.rep), pvs[, q], col="lightgrey")
  }
}

tval.thrsh.test = function(sample, qrange, theta0, sig=0.05, k_min=10){
  #' T metric independence test for threshold selection
  #'
  #' uses difference metric idea to check independence of X-u and Z, thus checking
  #' the suitability of qu as a threshold quantile
  #' does this for a range of quantiles, using an adaptive number of bands at each
  #' stage
  #'
  #' @param sample matrix of two columns containing sample data
  #' @param qrange candidate threshold quantiles for conditioning variable
  #' @param theta0 vector of initial H&T parameter values to be used in optimisation
  #' @param sig float level of sigficance to test against
  #' @param k_min int min number of points acceptable in each band
  #'
  #' @returns no value is returned
  #' @export
  pvalues=c()
  X = sample[,1] ; Y=sample[,1]

  N = length(X)
  k_seq = c()
  nbs = c()
  for (i_q in 1:length(qrange)){
    q = qrange[i_q]
    print(q)

    Ne = (1-q)*N
    n_bands_max = 40
    n_bands = min(n_bands_max,floor(Ne/k_min))

    k = Ne%/%n_bands
    k_seq[i_q] = k
    nbs[i_q] = n_bands

    print(n_bands)
    cond_test = ht_thrsh_vldtn(sample, q, theta0, plots = F, sig = sig,n_bands=n_bands, k_min=k_min)
    pvalues[i_q] = cond_test$p_value
  }

  plot(qrange, k_seq, pch=4, xlab="number in each band", ylab="quantile")
  lines(qrange, k_seq, col="blue")

  plot(qrange, pvalues, pch=4)
  abline(h=0.05, col="red")
  for (i_q in 1:length(qrange)){
    text(qrange[i_q], pvalues[i_q]+0.02, paste(nbs[i_q]))
  }

}

indv.thrsh.test = function(sample, qu, theta0, n_bands=5, sig=0.05, plots=T, k_min = 5, nr=1000){
  #' T value indepedence test for individual threshold
  #'
  #' uses difference metric idea to check independence of X-u and Z, thus checking
  #' the suitability of qu as a threshold quantile
  #'
  #' @param sample matrix of two columns containing sample data
  #' @param qu candidate threshold quantile for conditioning variable
  #' @param theta0 vector of initial H&T parameter values to be used in optimisation
  #' @param n_bands integer number of bands to split data into
  #' @param sig float level of sigficance to test against
  #' @param plots bool, true to show diagnostic plots
  #' @param k_min int min number of points acceptable in each band
  #' @param nr int number of randomised bands
  #'
  #'@returns list of resulting p-value boolean decision of accept/reject
  #'@export
  XandZ = get_X_and_Z(sample, qu, theta0)

  res_sample = matrix(nrow=length(XandZ$X), ncol=2)
  res_sample[,1] = XandZ$X ; res_sample[,2] = XandZ$Z

  bands = get_bands(res_sample, n_bands)

  if (length(bands$bandsY[[1]])<k_min){
    warning("too few elements in each band")
  }

  all_z = XandZ$Z
  cdfs = get_ecdfs(bands$bandsY)
  t_obs = get_t(cdfs)

  t_set = get_t_set(all_z, n_bands, nr)
  alpha = sig  # significance level of hypothesis test
  t_sig = quantile(t_set, 1- alpha)

  if(plots){
    par(mfrow=c(3,1))
    plot(XandZ$X, XandZ$Z, pch=16, xlab="X-u", ylab="Z", ylim=c(-10, 10))
    for (i in 1:n_bands){
      points(bands$bandsX[[i]], bands$bandsY[[i]], col=i+10, pch=16, cex=0.5)
      abline(v=max(bands$bandsX[[i]]), lwd=4, col="purple")
    }

    plot(ecdf(bands$bandsY[[1]]),main="", cex=0.5, xlim=c(-10,10), xlab="Z", ylab="empirical cdf", col=11)
    for (i in 2:n_bands){
      lines(ecdf(bands$bandsY[[i]]),main="",col=i+10, cex=0.5, xlab="Z", ylab="empirical cdf")
    }

    hist(t_set, xlim=c(min(t_set), max(max(t_set), t_obs)))
    abline(v=t_sig, col="red", lwd=2)
    abline(v=t_obs, col="blue", lwd=2)
    par(mfrow=c(1,1))
  }

  decision = as.numeric((t_sig[1] > t_obs))
  p_value = sum(t_set>t_obs)/nr

  return(list(p_value=p_value, decision=decision))

}


# internal functions ------------------------------------------------------


get_X_and_Z = function(sample, qu, theta0){
  #' get sorted X-u and Z for a given sample & fit conditions
  #'@keywords internal
  sample_x = sample[,1] ; sample_y = sample[,2]

  qu_fit = ht.fit(sample_x, sample_y, qu, keef=T, theta0=theta0, plot=F)

  u = quantile(sample_x, qu)
  excess_x = sample_x[sample_x>=u] - u

  s_excess_x = excess_x
  s_excess_x = sort(excess_x)
  s_res = c()

  for (i in 1:length(s_excess_x)){
    s_res[i] = qu_fit$res[excess_x==s_excess_x[i]]
  }

  return(list(X=s_excess_x, Z=s_res))

}

get_bands = function(sample, n_bands){
  #' split Y into given number of bands by X values
  #'@keywords internal
  X = sample[,1] ; Y = sample[,2]

  bands_X = eq_split(X, n_bands)$chunks
  bands_Y = eq_split(Y, n_bands)$chunks

  return(list(bandsX = bands_X, bandsY=bands_Y))

}

eq_split = function(X, n){
  #' splits a vector into n equal chunks (with overlap if necessary)
  #'@keywords internal
  k = length(X) %/% n
  X_overflow = X[(k*n+1): length(X)]

  chunks = list()
  for (i in 1:n){
    chunks[[i]] = X[(1+ (i-1)*k): (i*k)]
  }
  return(list(chunks=chunks, overflow=X_overflow))
}

get_ecdfs = function(bands){
  #' get ecdfs given bands
  #'@keywords internal
  n_bands = length(bands)
  
  emp_cdfs = list(knots(ecdf(bands[[1]])))

  for (i in 2:n_bands){
    emp_cdfs = list.append(emp_cdfs, knots(ecdf(bands[[i]])))
  }

  return(emp_cdfs)

}

get_t = function(cdfs){
  #' get test statistic t from quantiles and all residuals
  #'@keywords internal

  max_diffs = c()
  n_bands = length(cdfs)

  for (i in 1:n_bands){
    diffs = c()
    for (j in 1:n_bands){
      diffs[j] = max(abs(cdfs[[i]] - cdfs[[j]]))
    }
    max_diffs[i] = max(diffs)

  }

  t = max(max_diffs)

  return(t)

}

get_t_set = function(all_z, n_bands, k){
  #'@keywords internal
  t = c()

  for (i in 1:k){

    all_z_i = sample(all_z, length(all_z), replace=F)
    bands_i = eq_split(all_z_i, n_bands)$chunks
    cdfs_i = get_ecdfs(bands_i)
    t_i = get_t(cdfs_i)

    t[i] = t_i

  }

  return(t)

}
