#' build_interactome
#'
#' @param directory_interactome directory of String interactome file. ignored if
#'      interactome is not NULL
#' @param scores_threshold thresholds vector, if NULL no filtering will be applied
#' @param tax_id taxonomy id of the specie
#' @param interactome interactome data.frame, default = NULL.
#' @param AnnNbyN vector to split the searching on String. If the function gives
#'     Timeout Error try to lower the value (e.g 999 or 1499). Default = 1999.
#'
#' @importFrom openxlsx read.xlsx
#' @importFrom stringr str_split
#' @importFrom stats na.omit
#' @importFrom utils read.delim
#'
#' @return List containing information about interactome
#' @export
#'
build_interactome <- function(directory_interactome, tax_id,
                              interactome = NULL,
                              scores_threshold = NULL,
                              AnnNbyN = 1999){

  if (is.null(interactome)){
    interactome <- read.delim(directory_interactome, sep = " ")
  }

  interactome.filtered <- filter_interactome(interactome,
                                             scores_threshold = scores_threshold)


  all.interacting.proteins <- unique(c(interactome.filtered$protein1,
                                       interactome.filtered$protein2))

  ANN.MAP <- annotation_and_mapping(all.interacting.proteins, tax_id, AnnNbyN)

  ppi.int <- get_ppi_interactions(interactome.filtered,ANN.MAP$annotation.df,ANN.MAP$mapped_id)

  return(ppi.int)
}

