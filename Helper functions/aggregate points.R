#Aggregated data
point_aggregation <- function(points, boundary) {

  require(sf)

  # Ensure same CRS
  if (st_crs(points) != st_crs(boundary)) {
    stop("CRS of points and boundary must match.")
  }

  # Compute centroids of polygons
  centroids <- st_centroid(boundary)

  # Identify points inside any polygon
  inside_list <- st_within(points, boundary)
  inside_idx  <- lengths(inside_list) > 0

  # Points inside stay assigned to their polygon
  inside_assign <- rep(NA_integer_, nrow(points))
  inside_assign[inside_idx] <- sapply(inside_list[inside_idx], `[`, 1)

  # Points outside get assigned to nearest centroid
  outside_idx <- !inside_idx

  if (any(outside_idx)) {
    nearest_centroid_idx <- st_nearest_feature(points[outside_idx, ], centroids)
    inside_assign[outside_idx] <- nearest_centroid_idx
  }

  # Count number of points per polygon
  point_counts <- tabulate(inside_assign, nbins = nrow(boundary))

  return(point_counts)
}
