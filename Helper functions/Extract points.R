extract_points <- function(pp, crs) {

  require("sf", "spatstat", "dplyr")

  # Convert to data frame
  p_patterns <- data.frame(pp)

  # Create geometry object
  geometry_obj <- st_as_sf(p_patterns, coords = c("x", "y"), crs = crs) |>
    dplyr::mutate(id = seq_len(nrow(p_patterns)))

  return(geometry_obj)
}
