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
    # presence only thinning prob linear predictor
    #
    # h_s = covariates for thinning probability
    # cc  = presence-only data model regression coefficients
    #
    logit(b[pixel]) <-  inprod(h_s[pixel,] , cc)
  }
  # The presence only data model.
  #
  # This part of the model just uses the
  #  what we have calculated above (lambda
  #  and b). The denominator of this likelihood
  #  is actually a scalar so we can calculate it
  #  outside of a for loop. Let's do that first.
  #
  # The presence_only data model denominator, which
  #  is the thinned poisson process across the
  #  whole region (divided by the total number of 
  #  data points because it has to be 
  #  evaluated for each data point).
  po_denominator <- inprod(lambda[1:npixel], b[1:npixel] ) / m
  #
  # Loop through each presence-only datapoint
  #  using Bernoulli one's trick. The numerator
  #  is just the thinned poisson process for
  #  the ith data point.
  for(po in 1:m){
    ones[po] ~ dbern(
      exp(
        log(lambda[po_pixel[po]]*b[po_pixel[po]]) -
          po_denominator
      ) / CONSTANT
    )
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
# Priors for presence-only data model
for(po in 1:npar_po){
  cc[po] ~ dlogis(0, 1)
}
# Priors for det/non-det data model
for(pa in 1:npar_pa){
  a[pa] ~ dlogis(0, 1)
}
# Derived parameter, the number of cells occupied
zsum <- sum(z)
}