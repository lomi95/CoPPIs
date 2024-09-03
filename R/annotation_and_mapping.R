#' annotation_and_mapping
#'
#' @param all.interacting.proteins vector of proteins that need to be annotated
#' @param AnnNbyN number to split 'all.interacting.proteins' to retrieve annotation
#'     on StringDB
#' @param tax_id taxonomy id of the species
#'
#' @importFrom rbioapi rba_string_annotations
#' @importFrom rbioapi rba_string_map_ids
#'
#' @return A list with mapped ids data.frame and the annotation data.frame
#'
#' @export
annotation_and_mapping <- function(all.interacting.proteins, tax_id, AnnNbyN){

  vect.ann <- seq(1,length(all.interacting.proteins),by = AnnNbyN)

  all.annotation <- list()
  all.mapping <- list()
  for (i in 1:(length(vect.ann)-1)){
    all.annotation[[i]] <- rba_string_annotations(
      all.interacting.proteins[vect.ann[i]:(vect.ann[i+1]-1)],tax_id,verbose = F, split_df = F)
    all.mapping[[i]] <- rba_string_map_ids(
      all.interacting.proteins[vect.ann[i]:(vect.ann[i+1]-1)],tax_id,verbose = F)
  }
  all.annotation[[i+1]] <- rba_string_annotations(
    all.interacting.proteins[vect.ann[i+1]:length(all.interacting.proteins)],tax_id, split_df = F)
  all.mapping[[i+1]] <- rba_string_map_ids(
    all.interacting.proteins[vect.ann[i+1]:length(all.interacting.proteins)],tax_id)


  map_id.merged <- Reduce(rbind,all.mapping)

  annotation.df0 <- Reduce(merge_annotations,all.annotation)
  annotation.df1 <- annotation.df0[annotation.df0$number_of_genes>1,]

  return(list(mapped_id = map_id.merged,
              annotation.df = annotation.df1))
}
