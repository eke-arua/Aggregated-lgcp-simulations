points_to_boundary <- function(points_sf, boundary_sf, count = TRUE) {

  require(sf)
  require(dplyr)

  # Handle input types
  if (inherits(points_sf, "sf")) {
    points_sf <- st_transform(points_sf, st_crs(boundary_sf))
  } else if (inherits(points_sf, "ppp")) {
    df <- data.frame(x = points_sf$x, y = points_sf$y)
    points_sf <- st_as_sf(df, coords = c("x", "y"), crs = st_crs(boundary_sf))
  } else {
    points_sf <- st_as_sf(
      as.data.frame(points_sf),
      coords = c("x", "y"),
      crs = st_crs(boundary_sf)
    )
  }

  # Add region ID
  boundary_sf <- boundary_sf |>
    mutate(.region_id = dplyr::row_number())

  # Join points to polygons
  joined <- st_join(points_sf, boundary_sf, left = FALSE)

  if (count) {
    summary <- joined |>
      st_drop_geometry() |>
      group_by(.region_id) |>
      summarise(n_points = n(), .groups = "drop")

    boundary_sf |>
      left_join(summary, by = ".region_id") |>
      mutate(n_points = ifelse(is.na(n_points), 0, n_points))

  } else {
    return(joined)
  }
}
