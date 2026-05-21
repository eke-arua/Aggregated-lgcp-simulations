agg_counts <- function(boundary, points) {
  require("sf")

  if(st_crs(points) != st_crs(boundary)){
    stop("The crs coordinates do not match")
  }

  counts <- lengths(st_intersects(boundary, points))

  return(counts)

}
