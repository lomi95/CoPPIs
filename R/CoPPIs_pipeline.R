#' CoPPIs_pipeline
#'
#' @param dataset The dataset ProteinxSamples
#' @param genes_id genes id vector
#' @param names_of_groups names of the groups
#' @param name_analysis name of the analysis to be saved
#' @param interactome interactome
#' @param categories.interactome interactome for each category
#' @param annotation.interactome annotation of protein interactions
#' @param p.filter threshold CoPPIs pvalue, default 0.001
#' @param score.filter threshold CoPPIs score, default 2
#' @param gCOR.groups.l if a precedent analyisis was run, it could be used to speed up the computational time
#' @param gCOR.categories if a precedent analyisis was run, it could be used to speed up the computational time
#' @param signifCorr threshold significance correlation, default 0.05
#' @param correctionCorr correction method for pvalue correlation, default BH
#' @param min_n.corr threshold minimum observation pair for correlation, default NULL
#' @param significance.coppi threshold significance for CoPPIs pvalue, default 0.05
#' @param perc.ned.nppi threshold percentage number edges number ppi, default 0.01
#' @param perc.sign threshold percentage number edges significants on total edges, default 0.3
#' @param correction.coppi method for correction pvalues CoPPI, default BH
#' @param corr.test correlation test, pearson or spearman, default Spearman
#' @param compute_weights if T transformation formula will be applied to the non significant correlation according to their pvalues, default T
#' @param min_edges threshold, minimum edges of biological terms, default 1
#' @param max_edges threshold, maximum edges of biological terms, default Inf
#' @param categories name of the biological categories on to which apply CoPPIs, default  c("Component","Process","RCTM","WikiPathways")
#' @param threshold_edge_cyto threshold of similarity edges between terms to be loaded on Cytoscape, default 0.6
#' @param Cytoscape if T the result will be loaded on Cytoscape, default T
#'
#'
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph edge.attributes
#' @importFrom igraph intersection
#' @importFrom igraph E
#' @importFrom igraph intersection
#' @importFrom igraph intersection
#'
#' @importFrom Hmisc rcorr
#' @importFrom stringr str_to_title
#'
#' @importFrom RCy3 cytoscapePing
#' @importFrom RCy3 closeSession
#' @importFrom RCy3 saveSession
#'
#'
#' @return List file with every information about CoPPIs output
#' @export
#'
CoPPIs_pipeline <- function(dataset,
                            genes_id,
                            names_of_groups,
                            name_analysis,
                            interactome = interactome.hs,
                            annotation.interactome = annotation.interactome.hs,
                            categories.interactome = categories.interactome.hs,
                            p.filter = 0.001,
                            score.filter = 2,
                            gCOR.groups.l = NULL,
                            gCOR.categories = NULL,
                            signifCorr = 0.05,
                            correctionCorr = "BH",
                            min_n.corr = NULL,
                            significance.coppi = 0.05,
                            perc.ned.nppi = 0.01,
                            perc.sign = 0.3,
                            correction.coppi = "BH",
                            corr.test = "spearman",
                            compute_weights = T,
                            min_edges = 1,
                            max_edges = Inf,
                            categories = c("Component","Process","RCTM","WikiPathways"),
                            threshold_edge_cyto = 0.6,
                            Cytoscape = T){

  names_of_groups <- substr(gsub(" ",".",names_of_groups),1,12)
  name_analysis <- gsub(".xlsx","",name_analysis)
  if (!is.null(min_n.corr) && compute_weights){
    message("Warning: there is a threshold for correlation pairs \n\t the compute_weights setting will be set to FALSE")
  }

  if (Cytoscape) {
    tryCatch(cytoscapePing(), error = function(err) {
      message("Open Cytoscape, or restart it")
    })
  }
  rownames(dataset) <- genes_id
  if (max_edges<min_edges){
    message("max edges is lower than min edges, we are switching them")
    m1 <- max_edges
    max_edges <- min_edges
    min_edges <- m1
  }

  gInter <- graph_from_edgelist(as.matrix(interactome[,3:4]),directed = F)
  start.time.annotation <- Sys.time()
  prot.annotation <- lapply(annotation.interactome,function(x){
    compute_annotation(genes_id,x)
  })


  time_employed <- format_time_difference(difftime(Sys.time(),
                                                   start.time.annotation,units = "s"))
  message("Annotation computed in ", time_employed,sep = " ")

  names(prot.annotation) <- str_to_title(names(prot.annotation))
  categories <- str_to_title(categories)
  categories.avaliable <- names(prot.annotation)
  categ.notfound <- setdiff(categories,categories.avaliable)
  if (length(categ.notfound)){

    if (length(categ.notfound) == length(categories)){
      message("Error:\n",paste(setdiff(categories,
                                       categories.avaliable), ", "),
              "were not found in the avaliable categories")
      message("The categories avaliable are:\n", paste(categories.avaliable, ", "))
      return(Categories)
    } else if (length(categ.notfound) == 1){
      message("Warning:\n",paste(setdiff(categories,
                                         categories.avaliable), " "),
              "was not found in the avaliable categories")
      message("The categories avaliable are:\n", paste(categories.avaliable, ", "))
    } else {
      message("Warning:\n",paste(setdiff(categories,
                                         categories.avaliable), ", "),
              "were not found in the avaliable categories")
      message("The categories avaliable are:\n",
              paste(setdiff(categories.avaliable,
                            categories), ", "))
      categories <- intersect(categories, categories.avaliable)
    }
  }

  dataset.t <- t(dataset)
  ## correlazione
  names(names_of_groups) <- names_of_groups
  list.groups <- lapply(names_of_groups, function(x){
    d <- dataset.t[grep(x,rownames(dataset.t), ignore.case = T),]
    if (nrow(d)){
      return(d)
    } else {
      message("No columns named ",x," were found in the dataset")
    }
  })
  list.groups <- list.groups[!sapply(list.groups, is.null)]

  if (length(list.groups) < 2){
    message("Error: The number of groups are less than 2,
            please check the colnames(dataset)
            or the input names_of_groups")
    return(Groups)
  }


  if (is.null(gCOR.categories)){
    corr.coppi <- Correlations_CoPPIs(list.groups,
                                      names_of_groups,
                                      categories.interactome,
                                      gCOR.groups.l,
                                      min_n.corr,
                                      signifCorr,
                                      correctionCorr,
                                      corr.test,
                                      gInter,
                                      compute_weights)
    gCOR.categories <- corr.coppi$gCOR.categories
    gCOR.groups.l   <- corr.coppi$gCOR.groups.l

  }


  names(categories) <- categories

  results <- lapply(categories, function(x){
    names(categories.interactome) <- str_to_title(names(categories.interactome))
    return(CoPPI(annotations = prot.annotation[[x]],
                 correlations = gCOR.categories[[x]],
                 perc.ned.nppi = perc.ned.nppi,
                 perc.sign = perc.sign,
                 significance.coppi = significance.coppi,
                 correction.coppi = correction.coppi,
                 min_edges = min_edges,
                 max_edges = max_edges,
                 names_of_groups = names_of_groups,
                 p.filter = p.filter,
                 score.filter = score.filter,
                 interactome = categories.interactome[[x]]))
  })

  dir.create(name_analysis)
  setwd(name_analysis)
  tryCatch({
    lapply(results,save_results, name_analysis)
    cytoscapePing()
    closeSession(save.before.closing = F)
    all_excels <- dir()[grep(".xlsx",dir())]
    names(all_excels) <- all_excels

    lapply(all_excels,Load2Cytoscape,
           filter_similarity =  threshold_edge_cyto)
    setCytoStyle()
    saveSession(paste0(name_analysis,"_",threshold_edge_cyto,".cys"))
    message("Output saved in Cytoscape")
  }, error = function(err){
    message("Warning in Cytoscape:")
    message(err)
    message("\n")
  })


  CoPPI_results <- list(prot.annotation = prot.annotation,
                        gCOR.groups.l = gCOR.groups.l,
                        gCOR.categories = gCOR.categories,
                        resultsCoPPI = results,
                        parameters = list(name_analysis = name_analysis,
                                          names_of_groups = names_of_groups,
                                          signifCorr = signifCorr,
                                          correctionCorr = correctionCorr,
                                          min_n.corr = min_n.corr,
                                          perc.ned.nppi = perc.ned.nppi,
                                          perc.sign = perc.sign,
                                          significance.coppi = significance.coppi,
                                          correction.coppi = correction.coppi,
                                          corr.test = corr.test,
                                          compute_weights = compute_weights,
                                          min_edges = min_edges,
                                          max_edges = max_edges,
                                          categories = categories,
                                          threshold_edge_cyto = threshold_edge_cyto))
  save(CoPPI_results, file = paste0(CoPPI_results$parameters$name_analysis,".RDa"))
  setwd("..")
  return(CoPPI_results)
}
