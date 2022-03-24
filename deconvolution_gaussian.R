library('rstan')
library('shinystan')
source('stellar_rot_profile.R')

# read data
data = read.table('spec_cut_Ha.dat')
wave = data$V1
flux = data$V2

# define line a vsini
line = 6562.8
vsini = 50.0

# create rotational profile function
g_func = get_g(wave, line, vsini)

# create the input for stan, two gaussians
input_stan = list(N = length(wave), M=3, wave=wave-line, flux=flux, conv=g_func, 
                  sigma=0.004)

# run stan, 5,000 iters with adapt_delta of 0.85, take a bit
model_stan = stan(file='deconvolution_gaussian_gc.stan', data=input_stan, iter=4000,
                  cores=4, control = list(adapt_delta=0.80,
                                          max_treedepth=14), init=0)
# extract the parameters of the model
params = extract(model_stan)
a = params$a
b = params$b
c = params$c
flux_sim = params$flux_sim
flux_sim_conv = params$flux_sim_conv
flux_sim_rng = params$flux_sim_rng

plot_convoluted_rng_profiles(wave, flux, flux_sim_rng, 100)
median_param = get_median_profile_stan(a, b, c)
plot_unconvolved_profiles(wave-line, median_param, flux_sim, 100)
plot_convolved_profiles(wave-line, flux, g_func, median_param, flux_sim_conv, 100)

pars=c('a','b','c')
print(model_stan, pars=pars)

pairs(model_stan, pars=pars)
