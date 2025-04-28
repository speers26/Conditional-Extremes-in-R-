#
# Script containing functions for estimating bivariate densities using H&T (2004)
#

ht.pred.cnstr <- function(v, theta, Z, n_pred, dx, Y=T, points=F,
                          plot=T, adjust=0){
  #' Predicting HT points within set xrange
  #'
  #' Prediction of points using the ht(2004) model within a set range of x values
  #'
  #' @param v prediction threshold
  #' @param theta ht model parameters
  #' @param Z ht model residuals
  #' @param n_pred number of predictions to make
  #' @param dx width of restricted x range
  #' @param Y true if doing Y|X, false if doing X|Y
  #' @param points true to return prediction points in a matrix
  #' @param plot true to plot predicted points
  #' @param adjust inflation factor for residual noise sd, defaults to 0
  #' @returns matrix of predicted points
  #' @export

  set.seed(1)

  alpha <- theta[1] ; beta <- theta[2]

  rand_u = runif(n_pred) * pexp(dx)
  rand_exp = qexp(rand_u)

  rand_yi <- v + rand_exp
  rand_z <- sample(Z, n_pred, replace=T)

  bw = density(Z)$bw
  rand_z = rand_z + rnorm(length(rand_z), mean=0, sd=adjust*bw)

  prdctd_ymnsi <- rand_yi * alpha + rand_yi^beta * rand_z

  if(plot){
    if(1-Y){abline(h=v,col="darkred",lty="dashed")
      points(prdctd_ymnsi,rand_yi,col="darkred",pch=16,cex=0.5)}
    else{abline(v=v,col="darkred",lty="dashed")
      points(rand_yi,prdctd_ymnsi,col="darkred",pch=16,cex=0.5)}
  }

  if(points){
    return(matrix(data=c(rand_yi, prdctd_ymnsi), ncol=2))
  }
}

grid.dens = function(data, q.marginal, q.cond, xlim, ylim, nx, ny,
                     log=T, adjust=0){
  #' Gridded density estimate using H&T for extremes
  #'
  #' @description
  #' Calculates a grid based density estimate for bivariate data, by using the
  #' empirical density for regions in the body of the data and H&T (2004) for
  #' extreme regions
  #'
  #' GPD tails are fitted to marginal empirical distributions to transform to
  #' laplace margins
  #' @param data matrix of data in original margins
  #' @param q.marginal threshold exceedance quantile to use for margins
  #' @param q.cond threshold exceedance quantile to use for conditional model
  #' @param xlim x limits
  #' @param ylim y limits
  #' @param nx number of division points for xgrid
  #' @param ny number of division points for ygrid
  #' @param log set to true to plot log of density
  #' @param adjust inflation factor for residual noise sd, defaults to 0
  #'
  #' @returns dataframe of densities and probabilities
  #' @export

  # fit and transform margins -----------------------------------------------
  x = data[,1] ; y = data[,2]
  x.u = as.numeric(quantile(x, q.marginal))
  x.fit = gpd.fit(x, x.u, show=F)

  y.u = as.numeric(quantile(y, q.marginal))
  # y.fit = gpd.fit(y[x>7], y.u, show=T)
  y.fit = gpd.fit(y, y.u, show=T)

  max_y = y.u - y.fit$mle[1]/y.fit$mle[2]
  print(max_y)

  # change to Laplace margins
  x.l = qlaplace(as.numeric(as.vector(lapply(x, pspliced, x=x, u=x.u, gpd_par=x.fit$mle))))
  y.l = qlaplace(as.numeric(as.vector(lapply(y, pspliced, x=y, u=y.u, gpd_par=y.fit$mle))))

  ## remove infs
  y.l = y.l[is.finite((x.l))]
  x.l = x.l[is.finite((x.l))]
  x.l = x.l[is.finite(y.l)]
  y.l = y.l[is.finite(y.l)]


  # form grid --------------------------------------------------------------

  # get lines
  x.u.cond = quantile(x, q.cond, names=F)

  x.min = min(xlim)
  x.max = max(xlim)
  dx = (x.max - x.min)/nx
  x.div = seq(x.min, x.max, dx)
  x.mid = min(x.div[x.div>x.u.cond])

  y.min = min(ylim)
  y.max = max(ylim)
  dy = (y.max-y.min)/ny
  y.div = seq(y.min, y.max, dy)

  plot(x, y, pch=16, cex=0.5, xlim=xlim, ylim=ylim)

  abline(v=x.div, col="lightgrey", lty="dashed")
  abline(v=x.mid, col="blue", lty="dashed")
  abline(h=y.div,col="lightgrey", lty="dashed")

  # get centre points -------------------------------------------------------

  x.points = seq(x.min + 0.5*dx, x.max - 0.5*dx, by=dx)
  y.points = seq(y.min + 0.5*dy, y.max - 0.5*dy, by=dy)
  x.points.lower = x.points[x.points<x.mid]
  x.points.upper = x.points[x.points>x.mid]

  xy.grid.upper = expand.grid(x.points.upper, y.points)
  xy.grid.lower = expand.grid(x.points.lower, y.points)
  xy.grid = expand.grid(x.points, y.points)

  points(xy.grid.lower, pch=4, col="red", cex=0.3)
  points(xy.grid.upper, pch=4, col="green", cex=0.3)

  # form laplace grid --------------------------------------------------------

  y.div.l = qlaplace(as.numeric(as.vector(lapply(y.div, pspliced, x=y, u=y.u, gpd_par=y.fit$mle))))
  x.div.l = qlaplace(as.numeric(as.vector(lapply(x.div, pspliced, x=x, u=x.u, gpd_par=x.fit$mle))))

  x.mid.l = qlaplace(as.numeric(as.vector(lapply(x.mid, pspliced, x=x, u=x.u, gpd_par=x.fit$mle))))

  # get laplace points ------------------------------------------------------

  y.points.l = qlaplace(as.numeric(as.vector(lapply(y.points, pspliced, x=y, u=y.u, gpd_par=y.fit$mle))))
  x.points.l.lower = qlaplace(as.numeric(as.vector(lapply(x.points.lower, pspliced, x=x, u=x.u, gpd_par=x.fit$mle))))
  x.points.l.upper = qlaplace(as.numeric(as.vector(lapply(x.points.upper, pspliced, x=x, u=x.u, gpd_par=x.fit$mle))))

  y.points.l[is.infinite(y.points.l) & y.points.l > 0] = 1e10
  y.points.l[is.infinite(y.points.l) & y.points.l < 0] = -1e10
  x.points.l.lower[is.infinite(x.points.l.lower) & x.points.l.lower > 0] = 1e10
  x.points.l.lower[is.infinite(x.points.l.lower) & x.points.l.lower < 0] = -1e10
  x.points.l.upper[is.infinite(x.points.l.upper) & x.points.l.upper > 0] = 1e10
  x.points.l.upper[is.infinite(x.points.l.upper) & x.points.l.upper < 0] = -1e10

  xy.grid.l.upper = expand.grid(x.points.l.upper, y.points.l)
  xy.grid.l.lower = expand.grid(x.points.l.lower, y.points.l)

  # fix infs in laplace grid ------------------------------------------------

  x.div.l[is.infinite(x.div.l) & x.div.l > 0] = 2e10
  x.div.l[is.infinite(x.div.l) & x.div.l < 0] = -2e10
  y.div.l[is.infinite(y.div.l) & y.div.l > 0] = 2e10
  y.div.l[is.infinite(y.div.l) & y.div.l < 0] = -2e10

  # plotting laplace grid ---------------------------------------------------
  plot(x.l, y.l, cex=0.5, xlim=c(min(x.div.l[x.div.l>-1e10])-2, max(x.div.l[x.div.l<1e10])+2),
       ylim=c(min(y.div.l[y.div.l>-1e10])-2, max(y.div.l[y.div.l<1e10])+2), pch=16)
  abline(v=x.div.l, col="lightgrey", lty="dashed")
  abline(h=y.div.l, col="lightgrey", lty="dashed")
  # points(xy.grid.l.lower, pch=4, col="red", cex=0.3)
  points(xy.grid.l.upper, pch=4, col="green", cex=0.3)
  abline(v=x.mid.l, col="blue", lty="dashed")

  # get lower ps ------------------------------------------------------------
  probs.lower = apply(as.matrix(xy.grid.lower), 1, emp_box_prob, data=data,
                      x_div=x.div, y_div = y.div, x_mid_l = x.mid.l, x_mid = x.mid)
  #probs.lower = probs_lower/sum(probs_lower)*(1-0.5*exp(-hs_mid_l)) # quick fix for inconsistencies

  # get upper ps ------------------------------------------------------------
  ### fit model ###
  theta0 = c(0,0,1)

  fit = ht.fit(x.l, y.l, q.cond, keef=T, theta0=theta0)

  # doing for whole grid
  probs.upper = apply(as.matrix(xy.grid.l.upper), 1, box_prob, fit = fit,
                      x_div = x.div.l, y_div = y.div.l, adjust = adjust)
  # probs.upper = probs.upper/sum(probs.upper)*0.5*exp(-x.mid.l) # quick fix for inconsistencies

  # plot all probabilities --------------------------------------------------
  probs.spliced = c(probs.lower, probs.upper)
  grid = rbind(as.matrix(xy.grid.lower), as.matrix(xy.grid.upper))
  probs.df = data.frame(
    x = as.matrix(grid)[,1],
    y = as.matrix(grid)[,2],
    p = probs.spliced,
    dens = probs.spliced/(dx*dy)
  )

  if(log){
    fig = ggplot(probs.df, aes(x, y)) + geom_raster(aes(fill=log(dens))) +
    scale_fill_gradientn(colours=c("white", "yellow", "orange", "red", "black"), na.value = "white")
    fig = fig +
    ylim(layer_scales(fig)$y$range$range[1], layer_scales(fig)$y$range$range[2]) +
    theme_classic()
    fig = fig +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.01, 0.08), expand = c(0, 0)) + ylab("STP") + xlab("Hs")
  }
  else{
    fig = ggplot(probs.df, aes(x, y)) + geom_raster(aes(fill=(dens))) +
    scale_fill_gradientn(colours=c("white", "yellow", "orange", "red", "black"), na.value = "white")
    fig = fig +
    ylim(layer_scales(fig)$y$range$range[1], layer_scales(fig)$y$range$range[2]) +
    theme_classic() +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0))
  }

  print(fig)

  return(probs.df)
}


# internal functions ------------------------------------------------------

box_prob = function(point, fit, x_div, y_div, adjust=0){
  #' Estimates the probability of a point being in a box using HT (2004)
  #' 
  #' Uses a fitted conditional extremes (HT, 2004) model to estimate the probability of 
  #' lying within a box with given centre point
  #' 
  #' @param point the centre point of the box
  #' @param fit the fitted conditional extremes model
  #' @param x_div the division points for x
  #' @param y_div the division points for y
  #' @param adjust inflation factor for residual noise kernel smoothing sd, defaults to 0. Making this positive can reduce 'streak' effects in extrapolation
  #' 
  #' @keywords internal
  x_l = as.numeric(point[1])
  y_l = as.numeric(point[2])

  b_xl = max(x_div[x_l > x_div])
  b_xu = min(x_div[x_l < x_div])

  b_yl = max(y_div[y_l > y_div])
  b_yu = min(y_div[y_l < y_div])

  dx = b_xu - b_xl
  v = b_xl
  m = 100000
  pred_points = ht.pred.cnstr(v, fit$theta, fit$res, n_pred=m, dx=dx, points=T,
                              plot=F, adjust=adjust)

  cond_p_est = (1/2*exp(-b_xl)-1/2*exp(-b_xu)) * sum((b_xl<pred_points[,1] & pred_points[,1]<=b_xu) * (b_yl<pred_points[,2] & pred_points[,2]<=b_yu)) / m

  return(cond_p_est)

}

emp_box_prob = function(point, data, x_div, y_div, x_mid_l, x_mid){
  #' Counts proportion of points in box
  #'
  #' Takes the centre point of box and locations for division points, then calculates
  #' the proportion of points in those boxes
  #' 
  #' @param point the centre point of the box
  #' @param data the data to be used
  #' @param x_div the division points for x
  #' @param y_div the division points for y
  #' @param x_mid_l the box midpoitns for x on Laplace margins
  #' @param x_mid the box midpoint for x on original margins
  #' 
  #' @keywords internal
  x = as.numeric(point[1])
  y = as.numeric(point[2])

  b_xl = max(x_div[x > x_div])
  b_xu = min(x_div[x < x_div])

  b_yl = max(y_div[y > y_div])
  b_yu = min(y_div[y < y_div])


  emp_p = sum((data[,1] <= b_xu) * (data[,1] > b_xl) * (data[,2] <= b_yu) * (data[,2]> b_yl))/sum(data[,1]<x_mid)
  v = x_mid_l
  return(emp_p * (1-0.5*exp(-v)))
}

