library(sf)
library(here)
library(ggplot2)
library(rts2)
library(stars)
library(terra)
library(scampr)
library(dplyr)
library(parallel)
library(qs2)
library(spatstat)

#set seed
set.seed(290392)

#load shapefile for the boundary
b1 <- st_read(here("Reduced Cape Town.shp"))[1] |>
  dplyr::summarise(geometry = st_union(geometry)) |>
  st_transform(crs = 22234)

#load the boundary with the regions defined
b2 <- st_read(here("Reduced Cape Town.shp"))[1] |>
  st_transform(crs = 22234)

# create a grid object with rts2
g1 <- grid$new(b1, 700)

#Population density data from WorldPop
pop_density_stars <- read_stars("Pop density for reduced Cape Town.tif") |>
  st_transform(crs = 22234)

#Distance to healthcare facility
dist_stars <- read_stars("Distance to healthcare facility reduced.tif"  )/10000
dist_stars <- dist_stars |>
  st_transform(crs = 22234)

#Convert covariates to sf_objects
cape_raster <- st_as_sf(pop_density_stars, as_points = FALSE, merge = FALSE)
colnames(cape_raster) <- c("pop","geometry")

dist_raster <- st_as_sf(dist_stars, as_points = FALSE, merge = FALSE)
colnames(dist_raster) <- c("dist","geometry")

#Add the population density which is my offset to the grid
g1$add_covariates(cov_data = cape_raster,zcols = "pop")

#Add the distance and scaled distance covariate to the grid
g1$add_covariates(cov_data = dist_raster,zcols = "dist")
g1$grid_data$dist_scale <- scale(g1$grid_data$dist, center = TRUE, scale = TRUE)[,1]

#save the grid object (to be used later for fitting the models)
qs_save(g1, "Grid object.qs2")

#Convert the grid objects to rasters for simulations
boundary_raster <- terra::vect(g1$grid_data)

#Create a raster template

#res arguments should match that of grid to avoid aggregation bias (and possible)
#spatial misalignment
raster_template <- rast(boundary_raster, res = 700)

pop_raster <- rasterize(boundary_raster, raster_template, field = "pop")
dist_raster <- rasterize(boundary_raster, raster_template, field = "dist")
dist_raster_scale <- rasterize(boundary_raster, raster_template, field = "dist_scale")

#save these rasters for later use
# writeRaster(pop_raster, "Pop density reduced 700 meters.tif")
# writeRaster(dist_raster, "Distance to HCF reduced 700 meters.tif")
# writeRaster(dist_raster_scale, "Distance to HCF scaled reduced 700 meters.tif")



#Simulations
#true parameter values for the simulations
parameter_space <- list(
  beta_0 = c(-18, -17.75, -17.5),
  beta_1 = c(-3, -2, -1),
  scale = c(10000, 15000, 20000),
  var = 1) |>
  expand.grid()

#Sample from the parameter space
source("Helper functions/parameter sampling.R")

parameters <- parameter_sample(parameter_space, size = 1200)

#Simulations
sim_list <- mclapply(seq_len(nrow(parameters)), function(x) {

  params <- unlist(parameters[x, ])  # Extract row as vector

  #Linear predictor
  eta <- log(pop_raster) + params[1] + (params[2] * dist_raster_scale)
  intensity <- exp(eta)

  # Convert to an image format
  intensity_df <- as.data.frame(intensity, xy = TRUE)
  intensity_im <- vec2im(intensity_df$pop, intensity_df$x, intensity_df$y)

  # Simulate the log-Gaussian Cox process
  matern_sim <- spatstat.random::rLGCP(
    model = "matern",
    mu = log(intensity_im),
    var = params[4],
    scale = params[3],
    nu = 0.5,
    saveLambda = FALSE
  )

  return(matern_sim)
}, mc.cores = 4)

parameters$n_events <- sapply(sim_list, function(X) X$n)
parameters$index <- seq(nrow(parameters))

#save the parameter values
qs_save(parameters, "Parameters.qs2")

#Extract the points
points <- lapply(sim_list, function(X){
  as.data.frame(X) |>
    st_as_sf(coords = c("x", "y"), crs = 22234)

})

#Save the points
#qs_save(points, "Simulated points.qs2")

#Aggregate the disease cases to each boundary
source("Helper functions/points_to_boundary.R")

sim_list_agg <- lapply(sim_list, function(X){

  count = points_to_boundary(X, boundary_sf = b2)

  return(count$n_points)
})

#save the spatially aggregated counts
qs_save(sim_list_agg, "Aggregated simulations.qs2")
