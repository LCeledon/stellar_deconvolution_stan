# stellar_deconvolution_stan
Scripts and documents to perform stellar deconvolution using stan

This works includes, for the moment, the R scripts need to run deconvolution, as well as the profiles (absorption, emission and p-cygni) to test the deconvolution and the python script for fitting input gaussian.
It also includes the stan file required to run stan

It requires R with libraries 'rstan' and 'shynistan'.

deconvolution_ha_profiles.R : the main script that call stan in R via rstan

stellar_rot_profile.R : includes functions to create the stellar rotation profile, and diagnosis plot once the model has been run

deconvolution_gaussian_gc.stan : stan code, where the model is defined and then read and run by rstan

fit_gaussian_profile.ipynb : simple python script to fit a m-gaussian profile

OUT.HALPHA_VT010-Abs : test profile (unconvolved) for an absorption line around H\alpha

OUT.HALPHA_VT010-Emi : test profile (unconvolved) for an emission line around H\alpha

OUT.HALPHA_VT010-Pcygni : test profile (unconvolved) for a P-Cygni profile line around H\alpha
