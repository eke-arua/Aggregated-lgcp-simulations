parameter_sample <- function(df, size) {
  # Ensure the sample size is at least the number of rows in df
  if (size < nrow(df)) {
    stop("The 'size' argument cannot be less than the number of rows in the dataframe")
  }

  base_sample <- df

  # Sample additional rows with replacement
  extra_sample <- df[sample(nrow(df), size - nrow(df), replace = TRUE), ]

  final_sample <- rbind(base_sample, extra_sample)

  return(final_sample)
}
