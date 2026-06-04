library(rts2)
library(glmmrBase)
library(sf)
library(qs2)
library(dplyr)
library(parallel)


#Load the parameters used for the simulations
parameters <- qs_read("Parameters.qs2")

#load boundary
boundary_ct <- st_read("Reduced Cape Town.shp") |>
  st_transform(crs = 22234)

#Load grid with the covariates
grid_object  <- qs_read("Grid object.qs2")

#Initialize a grid for the aggregated models for rts2
grid_agg <- grid$new(boundary_ct, 700)

grid_agg$grid_data$dist_scale <- grid_object$grid_data$dist_scale
grid_agg$grid_data$pop <- grid_object$grid_data$pop

#Load simulated values
sim_list_agg <- qs_read("Aggregated simulations.qs2")

#fit lgcp models using rts2
results_agg <- mclapply(seq_along(sim_list_agg), function(i){

  #Add the simulated region counts to the grid object
  y_r <- sim_list_agg[[i]]

  grid_temp <- grid_agg$clone(deep = TRUE)
  grid_temp$region_data$y <- y_r

  print(paste0("running simulation ", i))

  fitg <- NULL

  #starting values to help fitting
  start_val <- if (parameters$scale[i] == 10000) {
    c(0, 9.2)
  } else if (parameters$scale[i] == 15000) {
    c(0, 9.6)
  } else if (parameters$scale[i] == 20000) {
    c(0, 9.9)
  } else {
    NULL
  }

  fitg <- tryCatch({
    if (is.null(start_val)) stop("Invalid scale")

    grid_temp$lgcp_ml(
      popdens     = "pop",
      covs        = "dist_scale",
      model       = "fexp",
      #iter_sampling = 1,
      # max_iter = 100,
      start_theta = start_val
    )

  }, error = function(e){
    message("Model failed (simulation ", i, "): ", e$message)
    return(NULL)
  })

  if (is.null(fitg) || is.null(fitg$coefficients) || nrow(fitg$coefficients) < 2) {
    return(
      matrix(NA_real_, nrow = 2, ncol = 7,
             dimnames = list(
               c("(Intercept)", "dist_scale"),
               c("par","est","SE","t","p","lower","upper")
             ))
    )
  }

  coefs <- fitg$coefficients[1:2, , drop = FALSE]
  return(coefs)
}, mc.cores = 6)

#Get models that have converged
keep <- !vapply(results_agg, function(x) all(is.na(x)), logical(1))

qs_save(keep, "Models that converged.qs2")
parameters_reduced <- parameters[keep, , drop = FALSE]

mod_results <- lapply(results_agg, as.data.frame)

mod_results <- lapply(mod_results, function(X){

  X$par = as.character(X$par)

  X
})

results_df <- dplyr::bind_rows(mod_results, .id = "ID")

#extract results for intercept and slope (fixed effects)
dist_res <- filter(results_df, par == "beta2")
intercept_res <- filter(results_df, par == "beta1")

#Function to help calculate coverage
in_between <- function(x, lower, upper) {
  x >= lower & x <= upper
}

#Function to convert the intercept (needed when spatstat is used to simulate)
b0_convert <- function(beta_0, grid_object){
  require(sf)

  beta_0 - log(as.numeric(st_area(grid_object$grid_data$bgrid[[1]])))
}



#Bias
bias_b0 <- mean((b0_convert(beta_0 = intercept_res$est, grid_object = grid_object) - parameters_reduced$beta_0))
bias_b1 <- mean((dist_res$est - parameters_reduced$beta_1))

#Relative Bias
rel_bias_b0 <- mean((b0_convert(beta_0 = intercept_res$est, grid_object = grid_object) - parameters_reduced$beta_0)/ parameters_reduced$beta_0)
rel_bias_b1 <- mean((dist_res$est - parameters_reduced$beta_1)/parameters_reduced$beta_1)

#Coverage
cov_b0 <- mean(in_between(parameters_reduced$beta_0, b0_convert(intercept_res$lower, grid_object), b0_convert(intercept_res$upper, grid_object)), na.rm = T)
cov_b1 <- mean(in_between(parameters_reduced$beta_1, dist_res$lower, dist_res$upper), na.rm = T)

#rsme
rmse_b0 <- sqrt(mean((b0_convert(intercept_res$est, grid_object) - parameters_reduced$beta_0)^2))


rmse_b1 <- sqrt(mean((dist_res$est - (parameters_reduced$beta_1))^2))

#montecarlo standard error
mcse_b0 <- sqrt((1/((nrow(parameters_reduced) - 1) * nrow(parameters_reduced))) * (sum((intercept_res$est - mean(intercept_res$est))^2)))
mcse_b1 <- sqrt((1/((nrow(parameters_reduced) - 1) * nrow(parameters_reduced))) * (sum((dist_res$est - mean(dist_res$est))^2)))

result_tab <- data.frame(rel_bias = c(rel_bias_b0, rel_bias_b1),
                         coverage = c(cov_b0, cov_b1),
                         rmse = c(rmse_b0, rmse_b1),
                         mcse = c(mcse_b0, mcse_b1),
                         row.names = c("Intercept", "Distance_scaled"))

#write.csv(result_tab, "results for aggregated point data simulations.csv")
write.csv(result_tab, "results for aggregated point data simulations aggregated_13_04_2026.csv")
