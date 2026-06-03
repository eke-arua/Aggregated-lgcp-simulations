library(sf)
library(exactextractr)
library(terra)
library(qs2)
library(spdep)

#covariates
dist_raster <- terra::rast("Distance to HCF reduced 700 meters.tif")
#dist_raster_scale <- terra::rast("Distance to HCF scaled reduced 700 meters.tif")
pop_raster <- rast("Pop density reduced 700 meters.tif")

#Load boundary
ct <- read_sf("Reduced Cape Town.shp")

#Create the covariates
ct$dist <- exact_extract(dist_raster, ct, 'mean')
ct$dist_scaled <- scale(ct$dist)[,1]
ct$pop <- exact_extract(pop_raster, ct, 'sum')


#Create a spatial weight matrix
nb <- poly2nb(ct)
W1 = nb2mat(nb, style = "B")


#Load the simulated aggregated realizations
sim_list <- qs_read("Aggregated simulations.qs2")

ct$Y <- sim_list[[1]]


S.CARbym(formula = Y ~ offset(log(pop)) + avg_dist,
       data = ct,
       family = "poisson",
       W = W1,
       burnin = 20000,
       n.sample = 50000,
       thin = 100)
