#' filter_interactome
#'
#' @param interactome data.frame or matrix interactome
#' @param scores_threshold a numeric vector named with the types of scores to filter.
#'     Default =  c("escore" = 0.15, "dscore" = 0.3)
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#'
#' @return A data.frame interactome filtered
#' @export
#'
filter_interactome <- function(interactome, scores_threshold = c("escore" = 150, "dscore" = 300)){

  score_type <- names(scores_threshold)

  if (sum(score_type %in% colnames(interactome))==length(score_type)){
    thres <- matrix(nrow = nrow(interactome),ncol = length(score_type))
    colnames(thres) <- score_type
    for (i in score_type){
      thres[interactome[,i] > scores_threshold[i],i] <- T
    }
    interactome <- interactome[rowSums(thres,na.rm = T)>=1,]
  } else {
    stop("names of 'score_threshold' are not matching with colnames(interactome)")
  }



  return(interactome)
}
