#' Title
#'
#' @param genes_id vector containing genes ids
#' @param annotation.df special data.frame containing annotation
#'
#' @return annotation data.frame with just and only the biological terms involved in the genes id list
#'
#'
compute_annotation <- function(genes_id,annotation.df){
  trova_geni <- function(x, genes_id) {
    genes_trovati <- intersect(x$preferredNames, genes_id)
  }

  genes_found <- apply(annotation.df,1,trova_geni,genes_id)
  annotation.df$genes_found <- genes_found
  annotation.df$number_of_genes_found <- sapply(genes_found, length)

  annotation.df.filt <- annotation.df[annotation.df$number_of_genes_found>1,]
  return(annotation.df.filt)
}
