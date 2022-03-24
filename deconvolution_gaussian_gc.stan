// Stan model to find the best M gaussian profile that convolved with the rotation functions
// best match the observed data
// para un futuro, stan vsini como parametro libre

functions{

 vector convolution(vector f, vector g, int N){
  vector[N] output;
  real suma;
  int condition;
  for (n in 1:N){
   suma = 0;
   for (k in 1:N){
    // integer division for N/2
    condition = n + k - 1 - N/2;
    if (condition > 1 && condition < N){
     suma += g[k] * f[condition];
    }
   }
   output[n] = suma;
  }
  return output;
 }
 
 real gauss(real x, real a, real b, real c){
  //just a gaussian profile
  return a * exp(-((x - b)/c)^2 / 2);
 }
}

data{
 int <lower=0> N; // len of input arrays
 int <lower=0> M; // number of gaussians to fit
 vector[N] flux; //input flux
 vector[N] wave; //input wavelength
 vector[N] conv; //input function that will convoluted later
 real <lower=0> sigma;
}

parameters{
 vector[M] a;
 ordered[M] b;
 vector<lower=0>[M] c;
}

transformed parameters{

 vector[N] flux_sim = rep_vector(0.0, N); //auxiliar vector that will store the gaussians
 vector[N] flux_sim_conv; //auxiliar vector convoluted with the conv profile
 
 for (n in 1:N){
  for (m in 1:M){
   flux_sim[n] += gauss(wave[n], a[m], b[m], c[m]);
  }
 }
 // the next line do the convolution
 flux_sim_conv = convolution(flux_sim, conv, N); 

}

model{

 //normal priors for a and b
 //gamma distribution for c
 for (m in 1:M){
  a[m] ~ normal(0.0,1.0);
  b[m] ~ normal(0.0,3.0);
  c[m] ~ gamma(2.0,1.0);
 }
  
 // the observed flux is distributed normally around the model
 flux ~ normal(flux_sim_conv, sigma);
}

generated quantities{
 real flux_sim_rng[N] = normal_rng(flux_sim_conv, sigma);
}
