#' Title
#'
#' @param x Number of signficative edges of a biological term
#' @param y Number of PPIs of the same biological term
#'
#' @return The CoPPIs score value
#' @export
#'
funscore <- function(x,y){
  R <- x/y
  if (R<=0.5){
    score <- log((1-R^2)/sqrt(R^3)*x^(5/2))
  } else {
    R1 <- 0.5
    score <- log((1-R1^2)/sqrt(R1^3)*(y/2)^(5/2)) +
      5/(2*x)*(x - y/2)
  }
  return(score)
}
