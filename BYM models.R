library(sf)
library(exactextractr)
library(terra)
library(qs2)
library(spdep)
library(CARBayes)
library(parallel)

#covariates
dist_raster <- terra::rast("Distance to HCF reduced 700 meters.tif")
#dist_raster_scale <- terra::rast("Distance to HCF scaled reduced 700 meters.tif")
pop_raster <- rast("Pop density reduced 700 meters.tif")

#Load boundary
ct <- read_sf("Reduced Cape Town.shp") |>
  st_transform(crs = 22234)

#Create the covariates
ct$dist <- exact_extract(dist_raster, ct, 'mean')
ct$dist_scaled <- scale(ct$dist)[,1]
ct$pop <- exact_extract(pop_raster, ct, 'sum')


#Create a spatial weight matrix
nb <- poly2nb(ct)
W1 = nb2mat(nb, style = "B")


#Load the simulated aggregated realizations
sim_list_agg <- qs_read("Aggregated simulations.qs2")

#Fit the BYM models
results <- mclapply(seq_along(sim_list_agg), function(i) {

  ct_i <- ct
  ct_i$Y <- sim_list_agg[[i]]

  mod.bym <- tryCatch(
    S.CARbym(
      formula = Y ~ offset(log(pop)) + dist_scaled,
      data = ct_i,
      family = "poisson",
      W = W1,
      burnin = 50000,
      n.sample = 550000,
      thin = 5000
    ),
    error = function(e) {
      message("Simulation ", i, " failed: ", conditionMessage(e))
      NULL
    }
  )

  return(mod.bym$summary.results[1:2, 1:3])


}, mc.cores = 12)

#Get the indices of the models that converged for the LGCP
keep <- qs_read("Models that converged.qs2")

#Extract results
mod_results <- results[keep]

mod_results <- lapply(mod_results, as.data.frame)

mod_results <- lapply(mod_results, function(X){

  X$par = ifelse(stringr::str_detect(row.names(X), "(Intercept)"), "Intercept", "dist_scaled")

  X
})

results_df <- dplyr::bind_rows(mod_results, .id = "ID")


#extract results for intercept and slope (fixed effects)
dist_res <- filter(results_df, par == "Intercept")
intercept_res <- filter(results_df, par == "dist_scaled")

#Function to help calculate coverage
in_between <- function(x, lower, upper) {
  x >= lower & x <= upper
}


#Bias
bias_b0 <- mean(intercept_res$Mean - parameters$beta_0)
bias_b1 <- mean(dist_res$Mean - parameters$beta_1)


#Coverage
cov_b0 <- mean(in_between(parameters$beta_0, intercept_res$`2.5%`,intercept_res$`97.5%`))
cov_b1 <- mean(in_between(parameters$beta_1, dist_res$`2.5%`, dist_res$`97.5%`))

#rsme
rmse_b0 <- sqrt(mean((intercept_res$Mean - parameters$beta_0)^2))
rmse_b1 <- sqrt(mean((dist_res$Mean - parameters$beta_1)^2))

#montecarlo standard error
mcse_b0 <- sqrt((1/((nrow(parameters) - 1) * nrow(parameters))) * (sum((intercept_res$Mean - mean(intercept_res$Mean))^2)))
mcse_b1 <- sqrt((1/((nrow(parameters) - 1) * nrow(parameters))) * (sum((dist_res$Mean - mean(dist_res$Mean))^2)))

result_tab <- data.frame(avg_bias = c(bias_b0, bias_b1),
                         #rel_bias = c(rel_bias_b0, rel_bias_b1),
                         coverage = c(cov_b0, cov_b1),
                         rmse = c(rmse_b0, rmse_b1),
                         mcse = c(mcse_b0, mcse_b1),
                         row.names = c("Intercept", "Distance_scaled"))

#write.csv(result_tab, "results for aggregated point data simulations.csv")
write.csv(result_tab, "results for BYM models_11_06_2026.csv")
