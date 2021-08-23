##############################################
#
# Simulate and analyze a poisson point process
#
# Written by M. Fidino
#
##############################################

# load the plotting and simulating functions
lapply(
  list.files(
    "./utilities/",
    full.names = TRUE
  ),
  source
)

# The packages we will need
packs <- c(
	"raster",
	"mvtnorm",
	"runjags", 
	"vioplot",
	"scales",
	"coda"
)

# this will load these packages (and download if necessary.)
package_load(packs)

set.seed(823)

# Whether or not you want to see the plots from each data
#  simulating step.
do_plots <- TRUE

# bounds of the plane we are sampling within
plane_xmin <- -1
plane_xmax <-  1
plane_ymin <- -1
plane_ymax <-  1

# number of pixels in space. 
npixels_x <- 1000
npixels_y <- 1000

# create a raster, currently has no data. 
#  We also make a blank raster
#  so that we can easily add other layers.
plane <- blank <- raster(
	ncol = npixels_x,
	nrow = npixels_y,
  xmn = plane_xmin,
	xmx = plane_xmax,
	ymn=plane_ymin,
	ymx=plane_ymax
)

# the x y coordinates of each pixel in the plane
plane_coord <- xyFromCell(
	plane,
	1:ncell(plane)
)

# generate a spatial covariate. 
#  gen_mvn from sim_utility.R script.
x_cov <- 0.4 * gen_mvn(plane, c(0.1,0.8), c(1,1), 0.1) + 
  0.6 * gen_mvn(plane, c(0.7,0.2), c(0.5,0.75), 0.4)

# add this covariate to plane raster object.
values(plane) <- as.numeric(scale(x_cov))
names(plane) <- 'x'

# Species presence across the landscape.
#  gen_process from ./utilities/sim_utility.R script.
# Inputs
#  rast = the raster with the spatial covariate
#  beta = matrix of regression coefficients, 
#    the intercept is the first value, and then
#    you need one for each covariate in dm.
#  my_seed = for set.seed() in function
#  dm = the design matrix
sp_pres <- gen_process(
	rast = plane,
	beta = matrix(
		c(6,-0.75),
		ncol = 2,
		nrow = 1,
		byrow=TRUE
	),
	my_seed = 4,
	dm = cbind(
		1, values(plane$x)
	)
)

# plot out species presence across landscape. 
#  plot_dist from ./utilities/plot_utility.R script
if(do_plots){
  plot_dist(
  	plane,
  	cov_name = "x",
  	pixel_id = sp_pres$pixel_id
  )
}

#############################
# simulate presence only data
#############################

# Calculating this at the same resolution as 
#  the latent Poisson point process.
#  We will later aggregate this to be 
#  at the same scale as the presence
#  absence data.

# make a detection covariate that influences presonce only detection
w_cov <- 0.5 * gen_mvn(plane,c(0.7,0.8), c(0.2,0.3), 0.2) +
  0.5 * gen_mvn(plane, c(0.1,0.3), c(0.2,0.3), 0.4)

# add this to the raster
temp <- blank
values(temp) <- as.numeric(scale(w_cov))
names(temp) <- "det"
plane <- addLayer(plane, temp)

# detection linear predictor for thinned Poisson process
beta_det <- c(1, 3)

# sample presence only from the true latent state
#
# rast = the study area
# dm = the presence only data model design matrix
# beta_det = regression coefficients for this model
# pres = where the species is present on the landscape
#  which gets sampled from for the po data.
po_data <- sample_po(
	rast = plane, 
	dm = cbind(1, values(plane$det)),
	beta_det,
	pres = sp_pres
)

# plot out just the presence only data
if(do_plots){
  plot_dist(
  	plane,
  	cov_name = "x",
  	pixel_id = po_data$pixel_id
  )
}

#################################
# Aggregate data to model it
#################################

# modify this based on size of your study
#  a larger agg factor will result in
#  fewer cells in the analysis 
agg_factor <- 30
agg_plane <- raster::aggregate(
	plane,
	fact = agg_factor
)

# function to aggregate pixels from output
#  of gen_process() and sample_po()
#  It aggregates the finer resolution
#  species presence data down to
#  the resolution we will model it at.
agg_po_pixel_id <- agg_pres(
	rast = plane,
	pixel_id = po_data$pixel_id, 
	agg_factor = agg_factor
)

# this is what the aggregated species 
#  presence in the region looks like
#  for the presence only data
if(do_plots){
  plot_dist(
  	agg_plane,
  	cov_name = "x",
  	pixel_id = agg_po_pixel_id
  )
}

##################################
# generate presence absence data
##################################

# We need to discretize the sample area to 
#  be more course (assuming that
#  our sampling method captures a slightly 
#  wider area than how we simulated
#  the latent Poisson process). 

agg_plane <- aggregate(
	plane,
	fact = agg_factor
)
agg_loc <- xyFromCell(
	agg_plane,
	1:ncell(agg_plane)
)

# We need to know which cells the species is 
#  and is not in this aggregated cell.
agg_pixelid <- agg_pres(
	plane,
	pixel_id = sp_pres$pixel_id,
	agg_factor
)

# sample the presence absence data
#
# rast = study area
# n = number of sites sampled
# visits = number of repeat surveys per site
# pixel_id = pixel locations of sites
# det_prob = detection probability given presence
pa_data <- sample_pa(
	rast = agg_plane, 
	n = 100,
	visits = 4, 
	pixel_id = agg_pixelid,
	det_prob = 0.3
)

if(do_plots){
  plot_dist(
  	agg_plane,
  	"x",
  	agg_pixelid
  )
  points(
  	agg_loc[pa_data$site_pixel,],
  	pch = 16,
  	col = "red"
  )
  par(xpd = NA)
  legend(
    x = 0.4, y = 1.3, pch = 16,
    col = c("black", "red"),
    legend = c("aggregated species presence", "site sampled"),
    bty = "n")
}

# fit an occupancy model with the grid based approach

# the number of grid points
G <- ncell(agg_plane)

# model using just the presence absence data
my_data <- list(
	# Total number of grid points
	npixel = G,
	# The occupancy covariates of 
	#  the latent state Poisson process
	x_s = cbind(1, values(agg_plane$x)),
	# The detection covariates of 
	#  the presence absence data
	v = matrix(1, ncol = 1, nrow = G),
	# The pixels the presence only data occurred
	pa_pixel = pa_data$y_mat$pixel,
	# The number of days a species 
	#   was detected per site
	#   for the presence absence data model
	y = pa_data$y_mat[,-1],
	# The number of presence absence sites
	nsite = nrow(pa_data$y_mat),
	# The log cell area.
	cell_area = log(prod(res(agg_plane))),
	# the number of latent parameters
	nlatent = ncol(sp_pres$beta),
	# The number of detection model parameters
	npar_pa = 1,
	# The number of repeat surveys
	W = 4
)

m1 <- run.jags(
	model = "./JAGS/occupancy_model.R", 
	data = my_data, 
	n.chains = 4, 
	inits = inits_pa, 
	monitor = c("beta", "a", "zsum"), 
	adapt = 1000, 
	burnin = 50000, 
	sample = 10000,
	thin = 2,
	method = 'parallel',
	modules = "glm",
	summarise = FALSE
)

# Fit the integrated model


my_data <- list(
	# Total number of grid points
	npixel = G, 
	# The occupancy covariates of the latent state Poisson process
	x_s = cbind(1, values(agg_plane$x)),
	# The detection covariates for the presence only data
	h_s = cbind(1, values(agg_plane$det)),
	# The detection covariates for the presence / absence data
	v = matrix(1, ncol = 1, nrow = G),
	# The pixels that the presence absence data occurred
	pa_pixel = pa_data$y_mat$pixel,
	# The pixels that the presence only data occurred							
	po_pixel = unlist(agg_po_pixel_id),
	# The number of presence only data points per season
	m = length(agg_po_pixel_id),
	# The number of days a species was detected per site / season
	#   for the presence absence data
	y = pa_data$y_mat[,-1],
	# The number of sites presence absence data was sampled
	nsite = nrow(pa_data$y_mat),
	# The log cell area, logged as the parameters are on the log scale
	#   in the model.
	cell_area = log(prod(res(agg_plane))),
	# Ones for the 'ones' trick in JAGS (for coding up likelihood)
	ones = rep(1, length(agg_po_pixel_id)),
	# A big constant value for 'ones' trick.
	CONSTANT = 10000,
	# Number of latent parameters
	nlatent = 2,
	# Number of observational parameters, presence only
	npar_po = 2,
	# Number of observational parameters, presence / absence
	npar_pa = 1,
	# Number of repeat surveys,
	W = 4
	)

m2 <- run.jags(
  model = "./JAGS/integrated_model.R", 
  data = my_data, 
  n.chains = 4, 
  inits = inits, 
  monitor = c("beta", "cc", "a", "zsum"), 
  adapt = 1000, 
  burnin = 50000, 
  sample = 20000,
  thin = 3,
  modules = "glm",
  method = 'parallel'
)


##########################
# summarize the two models
##########################

# get mcmc matrix from PA model
#  Removing the first column because it just signifies which chain the mcmc
#  step came from.
pamm <- as.matrix(as.mcmc.list(m1), chains = TRUE)[,-1]

# calculate median and 95% CI
pamm_ci <- apply(pamm, 2, quantile, probs = c(0.025,0.5, 0.975))

# do the same for the integrated model
intmm <- as.matrix(as.mcmc.list(m2), chains = TRUE)[,-1]
intmm_ci <- apply(intmm, 2, quantile, probs = c(0.025,0.5, 0.975))
windows(4,4)
tiff("model_comparison.tiff", height = 4, width = 4, units = "in",
		 res = 300, compression = "lzw")
par(mar = c(2,4,1,1))
plot(1~1, bty = "l", ylim = c(0, 20), xlim = c(0.5,2.5), type = "n",
		 xlab = "", ylab = "", xaxt = "n", yaxt = "n")

axis(1, at = seq(1, 2, 1), labels = F, tck = -0.025)
mtext(c("Intercept", "Slope"), 1, line = 0.6, at = c(1,2), cex = 1.5)

axis(2, at = seq(0,20,1), labels = F, tck = -0.025)
axis(2, at = seq(0,20,1/2), labels = F, tck = -0.025/2)
mtext(text = seq(0,20,1),2, line = 0.75, las = 1, at = seq(0,7,1))
mtext("Estimate", 2, line = 2, at = mean(c(0,7)), cex = 1.5)

lines(x = c(0,1.5), y = c(6,6), lty = 3, lwd = 3)
my_vioplot(pamm[,1], add = TRUE, wex = 0.3, at = 0.75, 
					 col = scales::alpha("#FF7E00", 0.4))
my_vioplot(intmm[,1], add = TRUE, wex = 0.3, at = 1.25, 
					 col = scales::alpha("#7e7F8B", 0.4))

lines(x = c(1.5,3), y = c(-0.75,-0.75), lty = 6, lwd = 3)

my_vioplot(pamm[,2], add = TRUE, wex = 0.3, at = 1.75, 
					 col = scales::alpha("#FF7E00", 0.4))
my_vioplot(intmm[,2], add = TRUE, wex = 0.3, at = 2.25, 
					 col = scales::alpha("#7e7F8B", 0.4))

legend("topright", legend = c("PA", "PA & PO"), 
			 col = c(scales::alpha("#FF7E00",0.4),scales::alpha("#7e7F8B",0.4)),
			 				lty = 1, lwd = 7 , bty = "n", title = "Model")
dev.off()



