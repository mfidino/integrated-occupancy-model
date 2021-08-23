##################################
#
# Utility functions for simulation
#
# Written by M. Fidino
# 
##################################


# this function loads pacakges and installs them if needed
package_load<-function(packages = NA, quiet=TRUE, verbose=FALSE, warn.conflicts=FALSE){
	
	# download required packages if they're not already
	pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
	if(length(pkgsToDownload)>0)
		install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", quiet=quiet, verbose=verbose)
	
	# then load them
	for(i in 1:length(packages))
		require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
}

# gen_mvn: 
#  simulate covariate values over cells in a raster
#
# rast = raster object of plane
#
# mu = numeric vector of length 2. Each mu is proportionally where you want
#  the mean value to be. (0,0) is the bottom left, (1,1) is top right.
# 
# sigma = Variance of covariate on x and y axes
#
# rho = Correlation of covariate on x and y axes
gen_mvn <- function(rast = NULL, mu = NULL,
										sigma = NULL, rho = NULL){
	# error checking
	if(length(mu) != 2 | !is.numeric(mu) | any(mu > 1) | any(mu < 0)){
		stop("mu must be a numeric vector of length 2.")
	}
	if(length(sigma) != 2 | !is.numeric(sigma) | any(sigma < 0)){
		stop("Sigma must be a non-negative numeric vector of length 2.")
	}
	if(length(rho) != 1 | !is.numeric(rho)| rho > 1 |rho < -1){
		stop("rho must be a numeric scalar between -1 and 1.")
	}
	
	# get bounds of raster
	bounds <- extent(rast)
	
	# input a proportion of where you want mu to be on x and y
	mu_loc <- c(bounds@xmin + mu[1] * (bounds@xmax - bounds@xmin),
							bounds@ymin + mu[2] * (bounds@ymax - bounds@ymin))
	
	Sigma <- diag(c(sigma[1] * abs(bounds@xmax - bounds@xmin),
									sigma[2] * abs(bounds@ymax - bounds@ymin)))
	# fill the off diagonal
	Sigma[2:3] <- rep(rho * prod(diag(Sigma)))
	
	to_return <- dmvnorm(xyFromCell(rast, 1:ncell(rast)), 
											 mean=mu_loc, 
											 sigma=Sigma)
}

# gen_process: 
#  Uses a PPP to place a species throughout a given landscape.
#  returns the given parameter values, seed used to generate distribution,
#  and the cell a species is located. Returns the pixels the species are located
#  as a vector for one sampling season or as a list if multiple sampling seasons.
#
# Arguments:
# 
# rast = raster object of spatial covariate
#
# dm = design matrix. First column must be vector of 1's for intercept. Must
#  have one row for each cell. Each column is for a different parameter.
#
# beta = conformable vector or matrix to values in dm object. If a matrix,
#   then each row is assumed to be a new point sampling session.
# 
# my_seed = seed to randomly generate values. For reproducibility.
#
# return_occ = return derived occupancy probability in each cell. 
#  defaults to false.

gen_process <- function(rast = NULL, dm = NULL, beta = NULL, my_seed = NULL,
												return_occ = FALSE){
	if(is.null(my_seed)){
		my_seed <- floor(runif(1, -1e4, 1e4))
	}
	set.seed(my_seed)
	
	# cell area on log scale
	cell_area <- log(prod(res(rast)))
	# coords
	rast_coord <- xyFromCell(rast, 1:ncell(rast))
	
	# log-linear intensity
	lin_pred <- t(exp(beta %*% t(dm)  + cell_area))
	
	# get area of raster
	rast_extent <- extent(rast)
	rast_area <- (rast_extent[2] - rast_extent[1]) * 
		(rast_extent[4] - rast_extent[3])
	
	# probability of an individual being in a cell
	my_prob <- 1 - exp(-lin_pred)
	
	# expected population size
	pop_size <- colSums(lin_pred)
	
	# add some variability to population size
	pop_size <- rpois(length(pop_size), pop_size)
	
	# place individuals on the landscape relative to the quality of each
	#  cell.
	pixel_id <- vector("list", length = length(pop_size))
	
	for(species in 1:length(pop_size)){
		pixel_id[[species]] <- sort(sample(x = 1:ncell(rast), 
																size = pop_size[species], 
																prob = my_prob[,species], 
																replace = FALSE))
		
	}
	
	# Drop the lists if only analyzing data for one year
	if(length(pop_size) == 1){
		pixel_id <- unlist(pixel_id)
	}
	to_return <- list(beta = beta,
										seed = my_seed,
										pixel_id = pixel_id,
										lambda = lin_pred)
	# return derived occupancy probability for each cell
	if(return_occ){
		to_return <- c(to_return, occ_prob)
	}
	return(to_return)
}


# sample pa:
#  

sample_pa <- function(rast = NULL, n = NULL, 
											det_prob = NULL, visits = NULL,
											pixel_id = NULL, my_seed = NULL){
	if(is.null(my_seed)){
		my_seed <- floor(runif(1, -1e4, 1e4))
	}
	set.seed(my_seed)
	# dimensions of raster
	w <- dim(rast)[1]
	h <- dim(rast)[2]
	
	# very rough and somewaht even spacing of cameras
	my_cams <- floor(seq(1, ncell(rast), length.out = n))
	
	# move half the cameras a little bit so they don't end up in a line
	jiggle_cam <- sort(sample(1:n, floor(n/2), replace = FALSE))
	
	my_cams[jiggle_cam] <- floor( my_cams[jiggle_cam] + median(diff(my_cams)/2))
	
	
	# keep all cams within realistic pixel range
	if(any(my_cams > ncell(rast))){
		my_cams[my_cams > ncell(rast)] <- ncell(rast)
	}
	
	if(!is.list(pixel_id)){
		pixel_id <- list(pixel_id)
	}
	ntime <- length(pixel_id)
	
	pixel_can_detect <- lapply(pixel_id, function(x){
		sort(x[which(x %in% my_cams)])
	})
	
	y_mat <- matrix(0, ncol = (1 + length(pixel_id)), nrow = n)
	y_mat[,1] <- sort(my_cams)
	colnames(y_mat) <- c("pixel", paste0("y", 1:length(pixel_id)))
	y_mat <- data.frame(y_mat)
	
	for(time in 1:ntime){
		tmp <- rbinom(length(pixel_can_detect[[time]]), size = visits, det_prob)
	y_mat[which(y_mat$pixel %in% pixel_can_detect[[time]]),(time+1) ] <- tmp
	}
	
	to_return <- list(y_mat = y_mat,
										visits = visits,
										site_pixel = my_cams,
										det_prob = det_prob,
										seed = my_seed)
	return(to_return)
}


# sample_po: sample presence only data
#  returns a vector of the pixels that a species is detected as 'present' 
#  based off of the 'det' covariate layer in object 'rast' and the 'beta' values.
#  Must also include the true presence of indviduals across the landscape that
#  was calculated via gen_process.
sample_po <- function(rast = NULL, dm = NULL, beta = NULL, 
											pres = NULL, my_seed = NULL, prop_lost = 0.2 ){
	if(is.null(my_seed)){
		my_seed <- floor(runif(1, -1e4, 1e4))
	}
	set.seed(my_seed)
	
	# get number of seasons sampled
	ntime <- ifelse(is.matrix(sp_pres$beta), nrow(sp_pres$beta), 1)
	
	# make det beta a matrix if we simulated data for more than one time step
	#  but only supplied a single group of covariate values
	#  i.e., sampled through time but things did not change in detection.
	if(ntime > 1 & !is.matrix(beta)){
		beta <- matrix(beta, ncol = length(beta), nrow = nrow(sp_pres$beta),
									 byrow = TRUE)
	}
	
	
	# cell area in log scale
	cell_area <- log(prod(res(rast)))
	# coords
	rast_coord <- xyFromCell(rast, 1:ncell(rast))
	
	#
	lp_det <- t(plogis(beta %*% t(dm)))
	# turn detection to 0 if species not in cell. 
	
	# for multiple seasons
	if(ntime > 1){
		for(time in 1:ntime){
			lp_det[-sp_pres$pixel_id[[time]],time] <- 0
		}
	} else { # for a single season
	lp_det[-sp_pres$pixel_id] <- 0
	}
	
	pixel_id <- rbinom(n = ncell(rast) * nrow(sp_pres$beta), 1, prob = lp_det)
	pixel_id <- matrix(pixel_id, nrow = nrow(lp_det), ncol = ncol(lp_det))
	pixel_id <- which(pixel_id > 0, arr.ind = TRUE)
	# make into a list
	pixel_id <- split(pixel_id[,1], factor(pixel_id[,2]))
	
	# drop 20% of the points
	pixel_id <- lapply(pixel_id, function(x){ sample(x, 
					 size = floor(length(x) * (1 - prop_lost)), replace = FALSE)})

	# remove list structure from pixel_id if one time step
	if(ntime == 1){
		pixel_id <- unlist(pixel_id)
	}
	
	to_return <- list(pixel_id = pixel_id, seed = my_seed, beta = beta)
	return(to_return)
}

# agg sp_pres
agg_pres <- function(rast = NULL, pixel_id = NULL, agg_factor = NULL){
	
	# make into list for t < 2
	if(!is.list(pixel_id)){
		pixel_id <- list(pixel_id)
	}
	ntime <- length(pixel_id)
	# reduce rast to one layer if layers > 1
	if(nlayers(rast) > 1){
		rast <- dropLayer(rast, 2:nlayers(rast))
	}
	temp <- rast
	
	# make new layers for the presence of the species
	#  in this raster layer
	for(time in 1:ntime){
	tmp_vals <- rep(0, ncell(temp))
	tmp_vals[pixel_id[[time]]] <- 1
	if(time == 1){ # overwrite first layer
		values(rast) <- tmp_vals
		names(rast) <- paste0("z", time)
	} else { # add the rest
		values(temp) <- tmp_vals
		names(temp) <- paste0("z",time)
		rast <- addLayer(rast, temp)
	}
	}
	
	# do aggregation. Aggregating for each column that we just generated
	#  agg_factor is supplied argument.
	agg_pres <- aggregate(rast,
												by = paste0("z", 1:ntime),
												fact = agg_factor, fun = sum)
	
	agg_loc <- xyFromCell(agg_plane, 1:ncell(agg_plane))

	agg_pixelid <- which(values(agg_pres)>0, arr.ind = TRUE)
	if(ntime == 1){
	  agg_pixelid <- matrix(
	    agg_pixelid,
	    nrow = length(agg_pixelid),
	    ncol = 2
	  )
	  colnames(agg_pixelid) <- c("row", "col")
	  agg_pixelid[,2] <- 1
	}
	agg_pixelid <- split(agg_pixelid[,1], factor(agg_pixelid[,2]))
	
	if(length(agg_pixelid) == 1){
		agg_pixelid <- unlist(agg_pixelid)
	}
	return(agg_pixelid)
}


# initial values for presence absence

inits_pa <- function(chain){
	gen_list <- function(chain = chain){
		list( 
			z = rep(1, G),
			beta = rnorm(2),
			a = rnorm(1),
			.RNG.name = switch(chain,
												 "1" = "base::Wichmann-Hill",
												 "2" = "base::Marsaglia-Multicarry",
												 "3" = "base::Super-Duper",
												 "4" = "base::Mersenne-Twister",
												 "5" = "base::Wichmann-Hill",
												 "6" = "base::Marsaglia-Multicarry",
												 "7" = "base::Super-Duper",
												 "8" = "base::Mersenne-Twister"),
			.RNG.seed = sample(1:1e+06, 1)
		)
	}
	return(switch(chain,           
								"1" = gen_list(chain),
								"2" = gen_list(chain),
								"3" = gen_list(chain),
								"4" = gen_list(chain),
								"5" = gen_list(chain),
								"6" = gen_list(chain),
								"7" = gen_list(chain),
								"8" = gen_list(chain)
	)
	)
}

# initial values for integrated sdm

inits <- function(chain){
	gen_list <- function(chain = chain){
		list( 
			z = rep(1, G),
			beta = rnorm(2),
			a = rnorm(1),
			cc = rnorm(2),
			.RNG.name = switch(chain,
												 "1" = "base::Wichmann-Hill",
												 "2" = "base::Marsaglia-Multicarry",
												 "3" = "base::Super-Duper",
												 "4" = "base::Mersenne-Twister",
												 "5" = "base::Wichmann-Hill",
												 "6" = "base::Marsaglia-Multicarry",
												 "7" = "base::Super-Duper",
												 "8" = "base::Mersenne-Twister"),
			.RNG.seed = sample(1:1e+06, 1)
		)
	}
	return(switch(chain,           
								"1" = gen_list(chain),
								"2" = gen_list(chain),
								"3" = gen_list(chain),
								"4" = gen_list(chain),
								"5" = gen_list(chain),
								"6" = gen_list(chain),
								"7" = gen_list(chain),
								"8" = gen_list(chain)
	)
	)
}



my_vioplot <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
												horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
												lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
												at, add = FALSE, wex = 1, drawRect = TRUE) 
{
	datas <- list(x, ...)
	n <- length(datas)
	if (missing(at)) 
		at <- 1:n
	upper <- vector(mode = "numeric", length = n)
	lower <- vector(mode = "numeric", length = n)
	q1 <- vector(mode = "numeric", length = n)
	q3 <- vector(mode = "numeric", length = n)
	med <- vector(mode = "numeric", length = n)
	base <- vector(mode = "list", length = n)
	height <- vector(mode = "list", length = n)
	baserange <- c(Inf, -Inf)
	args <- list(display = "none")
	if (!(is.null(h))) 
		args <- c(args, h = h)
	for (i in 1:n) {
		data <- datas[[i]]
		data.min <- min(data)
		data.max <- max(data)
		q1[i] <- quantile(data, 0.25)
		q3[i] <- quantile(data, 0.75)
		med[i] <- median(data)
		iqd <- q3[i] - q1[i]
		upper[i] <- min(q3[i] + range * iqd, data.max)
		lower[i] <- max(q1[i] - range * iqd, data.min)
		est.xlim <- c(min(lower[i], data.min), max(upper[i], 
																							 data.max))
		smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
																		 args))
		hscale <- 0.4/max(smout$estimate) * wex
		base[[i]] <- smout$eval.points
		height[[i]] <- smout$estimate * hscale
		t <- range(base[[i]])
		baserange[1] <- min(baserange[1], t[1])
		baserange[2] <- max(baserange[2], t[2])
	}
	if (!add) {
		xlim <- if (n == 1) 
			at + c(-0.5, 0.5)
		else range(at) + min(diff(at))/2 * c(-1, 1)
		if (is.null(ylim)) {
			ylim <- baserange
		}
	}
	if (is.null(names)) {
		label <- 1:n
	}
	else {
		label <- names
	}
	boxwidth <- 0.05 * wex
	if (!add) 
		plot.new()
	if (!horizontal) {
		if (!add) {
			plot.window(xlim = xlim, ylim = ylim)
			axis(2)
			axis(1, at = at, label = label)
		}
		# box()
		for (i in 1:n) {
			polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
							c(base[[i]], rev(base[[i]])), col = col, border = border, 
							lty = lty, lwd = lwd)
			if (drawRect) {
				lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
							lty = lty)
				rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
						 q3[i], col = rectCol)
				points(at[i], med[i], pch = pchMed, col = colMed)
			}
		}
	}
	else {
		if (!add) {
			plot.window(xlim = ylim, ylim = xlim)
			axis(1)
			axis(2, at = at, label = label)
		}
		# box()
		for (i in 1:n) {
			polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
																							rev(at[i] + height[[i]])), col = col, border = border, 
							lty = lty, lwd = lwd)
			if (drawRect) {
				lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
							lty = lty)
				rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
						 	boxwidth/2, col = rectCol)
				points(med[i], at[i], pch = pchMed, col = colMed)
			}
		}
	}
	invisible(list(upper = upper, lower = lower, median = med, 
								 q1 = q1, q3 = q3))
}

