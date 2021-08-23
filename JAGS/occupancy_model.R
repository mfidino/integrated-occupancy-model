model{
  # Bayesian version of the Koshkina (2017) model.
  #
  # The latent-state model
  for(pixel in 1:npixel){
    # latent state linear predictor
    #
    # x_s  = covariates for latent state
    # beta = latent state model regression coefficients
    # cell_area = log area of grid cell 
    #
    log(lambda[pixel]) <-inprod(x_s[pixel,], beta) + cell_area
    # Species presence in a gridcell as a Bernoulli trial
    z[pixel] ~ dbern(1 - exp(-lambda[pixel]))
  }
#
# Detection / non-detection data model
for(site in 1:nsite){
  # detection/non-detection data model linear predictor
  #
  #  v = detection covariates for the entire region B
  #  a = det/non-det data model regression coefficients
  #
  logit(rho[site]) <-inprod(v[pa_pixel[site], ],a)
  # The number of detections for site is a binomial
  #  process with Pr(rho[site]) | z = 1 with
  #  W sampling occasions per site.
  y[site] ~ dbin(
    z[pa_pixel[site]] * rho[site],
    W
  )
}
#
# Priors for latent state model
for(latent in 1:nlatent){
  beta[latent] ~ dnorm(0, 0.01)
}
# Priors for det/non-det data model
for(pa in 1:npar_pa){
  a[pa] ~ dlogis(0, 1)
}
# Derived parameter, the number of cells occupied
zsum <- sum(z)
}