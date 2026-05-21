library(here)
library(sf)
library(terra)
library(spatstat)
library(dplyr)


g1 <- grid$new(cape_raster, cellsize = 1200)

g1$points_to_grid(point_data = points)


pop_density_stars <- read_stars("Cape Town population density.tif")

park_dist_stars <- read_stars("Distance to healthcare facility.tif")

cape_raster <- st_as_sf(pop_density_stars, as_points = TRUE, merge = FALSE) %>%
  st_transform(crs = 22234)

dist_raster <- st_as_sf(park_dist_stars, as_points = TRUE, merge = FALSE) %>%
  st_transform(crs = 22234)

cov1 <- st_join(g1$grid_data, cape_raster)

pop_density_cov <- cov1 |>
  group_by(grid_id) |>
  summarise(total_pop = sum(`Cape Town population density.tif`, na.rm = TRUE))

cov2 <- st_join(g1$grid_data, dist_raster)

park_dist_cov <- cov2 |>
  group_by(grid_id) |>
  summarise(dist = mean(`Distance to healthcare facility.tif`, na.rm = TRUE))

covariates <- grid$new(cape_raster, 1200)

covariates$grid_data$pop_cov <- pop_density_cov$total_pop
covariates$grid_data$dist_cov <- park_dist_cov$dist


g1$add_covariates(covariates$grid_data,
                  zcols = "dist_cov",
                  popdens = "pop_cov",
                  weight_type = "area")


g1$reorder("minimax")

g1$lgcp_ml(popdens = "pop_cov",
           covs = "dist_cov",
           approx = "nngp",
           m = 15,
           iter_warmup = 150,
           iter_sampling = 50)

ct <- summarise(cape_town, boundary = st_union(geometry))

tmap_mode("plot")
tm_shape(ct_check) + tm_polygons() +
  tm_shape(ct) + tm_borders()

plot(g1$grid_data[4])


ct_check <- g1$grid_data[(is.nan(g1$grid_data$dist_cov)),]
