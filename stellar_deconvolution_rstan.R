# load libraries
library('rstan')
library('shinystan')
source('stellar_rot_profile.R')

# define parameters
line = 6562.8
vsini = 50
noise = 1e-3
m = 4
cores = 4
chains = 4
iterstan = 5000
warmupstan = 2000

# load the data
input_spec = read.table('../profiles/OUT.HALPHA_VT010-Abs')
input_wav = input_spec$V3
input_flux = input_spec$V5

# prepare data
length_spec = length(input_wav)
x = seq(from=input_wav[1], to=input_wav[length_spec],
        length.out=length_spec) - line
flux_interp = approx(input_wav-line, input_flux, xout=x)
y = flux_interp$y - 1

# create g func and convolved the y profile
g = get_g(x, line, vsini)
y_conv = get_convolution(y, g) + rnorm(length_spec, mean=0.0, sd=noise)

# export data so it can be fitted with python later
# remember to adapt the savename properly
savename = 'savename'
write.table(data.frame(x,y_conv), savename, row.names=FALSE, col.names=FALSE)

# create the initial values for the chains
# it uses the solution of python
# one listchain per chain
sd = 5e-3
listchain1=list(
  a=c() + rnorm(m, mean=0.0, sd=sd),
  b=c() + rnorm(m, mean=0.0, sd=sd),
  c=c() + rnorm(m, mean=0.0, sd=sd)
)
listchain2=list(
  a=c() + rnorm(m, mean=0.0, sd=sd),
  b=c() + rnorm(m, mean=0.0, sd=sd),
  c=c() + rnorm(m, mean=0.0, sd=sd)
)
listchain3=list(
  a=c() + rnorm(m, mean=0.0, sd=sd),
  b=c() + rnorm(m, mean=0.0, sd=sd),
  c=c() + rnorm(m, mean=0.0, sd=sd)
)
listchain4=list(
  a=c() + rnorm(m, mean=0.0, sd=sd),
  b=c() + rnorm(m, mean=0.0, sd=sd),
  c=c() + rnorm(m, mean=0.0, sd=sd)
)

# create stan input
initstan = list(listchain1, listchain2, listchain3, listchain4)
inputstan = list(N=length(x), M=m, wave=x, flux=y_conv, conv=g, sigma=noise)
controlstan = list(adapt_delta=0.80, max_treedepth=15)

# run stan
model_stan = stan(file='deconvolution_gaussian_gc.stan', data=inputstan, 
                  iter=iterstan, warmup=warmupstan, chains=chains, cores=cores,
                  control=controlstan, init=initstan)

#save/load model if desired
savename = 'savename'
loadname = 'finished_models/abs_4g_vsini_200_noise_1e-2'
save(model_stan, file=savename)
load(file=loadname)

# print stan results
pars=c('a','b','c')
print(model_stan, pars=pars)
pairs(model_stan, pars=pars)

# extract parameters
params=extract(model_stan)
a=params$a
b=params$b
c=params$c
flux_sim=params$flux_sim
flux_sim_conv=params$flux_sim_conv
flux_sim_rng=params$flux_sim_rng

# export parameters if desired
paramsname = 'abs_4g_vsini_200_noise_1e-2_params.csv'
write.csv(params, file=paramsname)

# get median parameters from a,b and c
median_param = get_median_profile_stan(a, b, c)

# diagnosis plots
plot_convoluted_rng_profiles(x, y_conv, flux_sim_rng, 1000)
plot_unconvolved_profiles(x, y, median_param, flux_sim, 1000)
plot_convolved_profiles(x, y_conv, g, median_param, flux_sim_conv, 1000)

# sqr difference sum
sqr_dif_sum = sqr_sum_stan(x, median_param, y)

# equivalent width difference
ew_dif = equivalent_width_dif(x, median_param, y)