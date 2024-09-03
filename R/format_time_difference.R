#' Title
#'
#' @param time_diff function that rewrites the time difference
#'
#' @return A string
#'
#'
format_time_difference <- function(time_diff) {
  if (time_diff < 1) {
    return(paste0(round(time_diff * 1000, 2), " ms"))
  } else if (time_diff < 60) {
    return(paste0(round(time_diff, 2), " secs"))
  } else if (time_diff < 3600) {
    return(paste0(round(time_diff / 60, 2), " mins"))
  } else {
    return(paste0(round(time_diff / 3600, 2), " hours"))
  }
}
