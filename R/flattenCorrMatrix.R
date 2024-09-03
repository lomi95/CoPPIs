#' flattenCorrMatrix
#'
#' @param cormat matrix of correlation values
#' @param pmat   matrix of correlation pvalues
#' @param nmat   matrix of observation pairs
#'
#' @return A data.frame that flattens the correlation matrices
#' @export
#'
flattenCorrMatrix <- function(cormat, pmat, nmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = cormat[ut],
    p = pmat[ut],
    n = nmat[ut]
  )
}
