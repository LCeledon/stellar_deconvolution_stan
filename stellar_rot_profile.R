

# Functions created from the code in Python by Lientur

G_function <- function(x, delta, epsilon){
  y = 1-(x/delta)^2    # transformation of wavelengths
  G = (2*(1-epsilon)*sqrt(y) + 0.5*pi*epsilon*y)/(pi*delta*(1-epsilon/3))
  # remove nans and set the total sum to 1
  G[is.nan(G)] = 0
  G = center_g(G)
  return (G / sum(G))
}


center_g <- function(g){
  # center the profile g, so that maximum is placed at the middle of the array
  n = length(g)
  new_g = replicate(n, 0)
  not_zero = g[g > 0]
  m = length(not_zero)
  ll = n %/% 2 - length(not_zero) %/%  2 - 1
  j = 1
  for (i in 1:n) {
    if(i > ll + 1 && j < m + 1){
      new_g[i] = not_zero[j]
      j = j + 1
    }
  }
  return (new_g)
}

get_g <- function(x, line, vsini, epsilon = 0.6) {
  delta = (line*vsini)/299792
  return (G_function(x, delta, epsilon))
}

m_gaussians <- function(x,a,b,c) {
  # create a profile of m gaussians, which are defined by the vectors a, b and c
  y = rep(0.0, times=length(x))
  M = length(a)
  for (m in 1:M){
    y = y + (a[m] * exp(-((x-b[m])/c[m])^2 / 2))
  }
  return (y)
}

get_generated_quantities <- function(sim_data){
  # returns the median values at each point of input x for generated quantities
  N = dim(sim_data)[2]
  output_array = rep(0.0, times=N)
  for (n in 1:N){
    x = sim_data[,n]
    output_array[n] = median(x)
  }
  return (output_array)
}

get_convolution <- function(x, g){
  # convolutes a profile with the g function
  N1 = length(x)
  N2 = length(g)
  N_h = N2%/%2
  y = convolve(x, g, type='open')
  return (y[(1+N_h):(N_h+N1)])
}

plot_generated_quantities <- function(x, y, y_sim){
  # do a plot of the profile and compare with the original
  plot(x, y, col='black', type='l')
  lines(x, y_sim, col='blue', type='l')
  legend('topleft', legend=c('Original','Simulated'), col=c('black','blue'),
         lty=1:1)
}

plot_gaussian_stan <- function(x, a, b, c, n_random){
  # do a plot of the stan output without convolution
  N = dim(a)[1]
  M = dim(a)[2]
  a_median = rep(0.0, times=M)
  b_median = rep(0.0, times=M)
  c_median = rep(0.0, times=M)
  for (m in 1:M){
    a_median[m] = median(a[,m])
    b_median[m] = median(b[,m])
    c_median[m] = median(c[,m])
  }
  x_big = seq(from=x[1], to=x[length(x)], length.out=100)
  y_big = rep(0.0, times=100)
  for (m in 1:M){
    aux = gaussian(x_big, a_median[m], b_median[m], c_median[m])
    y_big = y_big + aux
  }
  plot(x_big, y_big, col='black', type='l', xlab='Wavelength', ylab='Flux - 1')
  # select number of random models and plot it too
  stopifnot(n_random > 0)
  args = 1:N
  args_random = sample(args, n_random)
  for (n in 1:n_random){
    n_ran = args_random[n]
    y_ran = rep(0.0, times=100)
    for (m in 1:M){
      y_ran = y_ran + gaussian(x_big, a[n_ran,m], b[n_ran,m], c[n_ran,m])
    }
    lines(x_big, y_ran, col='blue', type='l')
  }
  legend('topleft', legend=c('Median'), col=c('black'), lty = 1)
}


plot_convoluted_rng_profiles <- function(wave, flux, flux_sim, n_ran){
  # first plot the input profile
  # and then random models if n_ran is greater than 0
  
  stopifnot(n_ran > 0)
  
  N_mod = dim(flux_sim)[1]
  args = 1:N_mod
  args_sample = sample(args, size=n_ran)
  
  plot(wave, flux, col='black', type='l', lwd=1, main='Convoluted Profile',
       xlab='Delta Wavelength', ylab='Normalized Flux - 1')
  for (n in 1:n_ran){
    arg_i = args_sample[n]
    y_i = flux_sim[arg_i,]
    lines(wave, y_i, col=rgb(red=0.0, green=0.0, blue=1.0,alpha=0.25), type='l', 
          lwd=1)
  }
  lines(wave, flux, col='black', type='l', lwd=5)
  legend('bottomleft', legend=c('Input Profile','Stan RNG Quantities'), col=c('black','blue'),
         lty=1)
}

get_median_profile_stan <- function(a, b, c){
  # compute the median values for a, b and c and returns the median profile
  M = dim(a)[2]
  output_m = matrix(nrow=M, ncol=3)
  for (m in 1:M){
    a_ = a[,m]
    output_m[m,1] = median(a_)
    b_ = b[,m]
    output_m[m,2] = median(b_)
    c_ = c[,m]
    output_m[m,3] = median(c_)
  }
  return (output_m)
}

plot_unconvolved_profiles <- function(x, y, params_matrix, flux_sim, n_random){
  # plot the gaussian profile for median parameters
  # plot random unconvolved profiles
  
  y_median = m_gaussians(x, params_matrix[,1], params_matrix[,2], params_matrix[,3])
  plot(x, y, col='blue', type='l', lwd=1, main='Absorption True Profile', 
       xlab='Delta Wavelength', ylab='Normalized Flux - 1')
  
  N = dim(flux_sim)[1]
  args = 1:N
  args_sample = sample(args, n_random)
  for (i in 1:n_random){
    arg_i = args_sample[i]
    y_i = flux_sim[arg_i,]
    lines(x, y_i, col=rgb(red=0, green=0, blue=1.0, alpha=0.25), type='l', lwd=1)
  }
  lines(x, y_median, col='green', type='l', lwd=3)
  lines(x, y, col='black', type='l', lwd=3)
  legend('bottomleft', legend=c('True Profile','Stan Median Profile','Stan Samples Profiles'),
         col=c('black','green','blue'), lty=1)
}

plot_convolved_profiles <- function(x, y, g_func, params_matrix, flux_sim_conv,
                                    n_random){
  # plot the input profile
  # plot the gaussian profile for median parameters
  # plot random unconvolved profiles
  
  y_median = y_median = m_gaussians(x, params_matrix[,1], params_matrix[,2], params_matrix[,3])
  y_median_conv = get_convolution(y_median, g_func)
  plot(x, y, col='black', type='l', lwd=1, main='Convoluted Profile',
       xlab='Delta Wavelength', ylab='Normalized Flux - 1')
  
  N = dim(flux_sim_conv)[1]
  args = 1:N
  args_sample = sample(args, n_random)
  for (i in 1:n_random){
    arg_i = args_sample[i]
    y_i = flux_sim_conv[arg_i,]
    lines(x, y_i, col=rgb(red=0, green=0, blue=1.0, alpha=0.25), type='l', lwd=1)
  }
  lines(x, y_median_conv, col='green', type='l', lwd=3)
  lines(x, y, col='black', type='l', lwd=3)
  legend('bottomleft', legend=c('Input Profile','Stan Median Profile','Stan Sample Profile'),
         col=c('black','green','blue'), lty=1)
}

get_mode_profile_stan <- function(a, b, c){
  # compute the mode of the parameters a,b and c
  # use the hist func with number of bins determined by the Freedman-Diaconis rule
  # to get the number of bins and passed to hist function
  mgauss = dim(a)[2]
  output_m = matrix(nrow=mgauss, ncol=3)
  for (m in 1:mgauss){
    a_ = a[,m]
    output_m[m,1] = get_mode_from_hist(a_)
    b_ = b[,m]
    output_m[m,2] = get_mode_from_hist(b_)
    c_ = c[,m]
    output_m[m,3] = get_mode_from_hist(c_)
  }
  return (output_m)
}

freedman_diaconis_rule <- function(arr){
  # given an array, use the Freedman-Diaconis rules to return the number of bins
  iqr = IQR(arr)
  bin_width = 2*iqr*(length(arr)**(-1/3))
  bin_number = (max(arr)-min(arr))/bin_width
  return(ceiling(bin_number))
}

get_argmax <- function(arr){
  narr = length(arr)
  argmax = 0
  maxval = -1e100
  for (n in 1:narr){
    if (arr[n] >= maxval){
      maxval = arr[n]
      argmax = n
    }
  }
  return(argmax)
}

get_mode_from_hist <- function(arr){
  nbins = freedman_diaconis_rule(arr)
  arrhist = hist(arr, breaks=nbins, plot=FALSE)
  histargmax = get_argmax(arrhist$counts)
  histmode = arrhist$mids[histargmax]
  return(histmode)
}

sqr_sum_stan <- function(x, median_param, y_true){
  # get the median_param profile and compare with true profile
  nlength = length(x)
  y_median = m_gaussians(x, median_param[,1], median_param[,2], median_param[,3])
  sqr_dif_sum = 0
  for (n in 1:nlength){
    aux = (y_true[n] - y_median[n])**2
    sqr_dif_sum = sqr_dif_sum + aux
  }
  return (sqr_dif_sum)
}

equivalent_width <- function(x, y){
  # compute the equivalent width of the line.
  # it assumes is already normalized at zero level
  # uses the composite trapezoidal rule to perform the integration
  nlength = length(x)-1
  int = 0
  for (n in 1:nlength){
    int = int + (x[n+1]-x[n])*(y[n+1]+y[n])
  }
  return (-0.5*int)
}

equivalent_width_dif <- function(x, median_param, y_true){
  y_median = m_gaussians(x, median_param[,1], median_param[,2], median_param[,3])
  ew1 = equivalent_width(x, y_true)
  ew2 = equivalent_width(x, y_median)
  return (ew1 - ew2)
}
