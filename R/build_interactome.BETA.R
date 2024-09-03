#' build_interactome.BETA
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
build_interactome.BETA <- function(directory_interactome, tax_id,
                              interactome = NULL,
                              scores_threshold = c("experimental" = 150, "database" = 300),
                              AnnNbyN = 1999){

  if (is.null(interactome)){
    interactome <- read.delim(directory_interactome, sep = " ")
  }

  interactome.filtered <- filter_interactome(interactome,
                                             scores_threshold = scores_threshold)


  all.interacting.proteins <- unique(c(interactome.filtered$protein1,
                                       interactome.filtered$protein2))

  ANN.MAP <- annotation_and_mapping(all.interacting.proteins, tax_id, AnnNbyN)


  Corum.downloaded <- read.xlsx("C:/Users/WKS/Desktop/Coppi_new/CORUM download 2022_09_12.xlsx")
  Corum.organisms <- split(Corum.downloaded,Corum.downloaded$Organism)
  Corum.human <- Corum.organisms$Human
  gene_name.complexes <- str_split(Corum.human$`subunits(Gene.name)`, ";")
  gene_name <- unique(unlist(str_split(Corum.human$`subunits(Gene.name)`, ";")))

  AnnMapCORUM <- annotation_and_mapping(gene_name,9606,AnnNbyN)

  mapGene.name1 <- AnnMapCORUM$mapped_id
  mapGene.name  <- mapGene.name1[na.omit(match(all.interacting.proteins,mapGene.name1$stringId)),]
  Corum.light <- data.frame(term = Corum.human$ComplexID,
                            category = rep("CORUM", length(Corum.human$ComplexID)),
                            description = Corum.human$ComplexName)
  inputGenes <- lapply(gene_name.complexes, function(x){
    mapGene.name$stringId[na.omit(match(x,mapGene.name$queryItem))]
  })


  preferredName <- lapply(gene_name.complexes, function(x){
    mapGene.name$preferredName[na.omit(match(x,mapGene.name$queryItem))]
  })
  Corum.light$inputGenes      <- inputGenes
  Corum.light$preferredNames  <- preferredName
  Corum.light$number_of_genes <- sapply(preferredName, length)

  Corum.light1 <- Corum.light[Corum.light$number_of_genes>1,]


  annotation.df2 <- rbind.data.frame(ANN.MAP$annotation.df,
                                     Corum.light1)

  map_id.merged <- rbind(ANN.MAP$mapped_id, mapGene.name)


  ppi.int <- get_ppi_interactions(interactome.filtered,annotation.df2,map_id.merged)

  return(ppi.int)
}

