#' Title
#'
#' @param Graph_pathways.groups graph pathways object
#' @param g1 group selected
#'
#' @return a data.frame
mean_abs_cor <- function(Graph_pathways.groups,g1){
  asd <- sapply(Graph_pathways.groups[[g1]],function(x){
    if (length(E(x)>0)){
      mean(abs(edge.attributes(x)$weight))
    } else {
      return(0)
    }
  })
}
