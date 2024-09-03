#' get_ppi_interactions
#'
#' @param interactome interactome object
#' @param annotation.df annotation dataframe
#' @param map_id.merged mapped id dataframe
#'
#' @importFrom igraph induced_subgraph
#' @return A information list for interactome
#'
#'
get_ppi_interactions <- function(interactome,
                                 annotation.df,
                                 map_id.merged){

  p1.mapped <- map_id.merged$preferredName[match(interactome$protein1,map_id.merged$queryItem)]
  p2.mapped <- map_id.merged$preferredName[match(interactome$protein2,map_id.merged$queryItem)]

  interactome.final <- cbind.data.frame(interactome[,1:2],
                                        p1.mapped,p2.mapped,
                                        interactome[,3:10])
  g.inter.final.hs <- graph_from_edgelist(as.matrix(interactome.final[,3:4]),
                                          directed = F)

  ann.splitted <- split(annotation.df,annotation.df$category)
  Gp.interactome <- lapply(ann.splitted,function(y){
    gp.x <- apply(y,1,function(x){
      tryCatch({
        return(induced_subgraph(g.inter.final.hs,x$preferredNames))
      }, error = function(e) {
        message("something went wrong, check your interactome\n")
        message(e)
        return(0)
      })
    })
  })

  for (i in 1:length(Gp.interactome)){
    names(Gp.interactome[[i]]) <- ann.splitted[[i]]$description
  }

  Proc.length <- lapply(Gp.interactome,function(x) {
    sapply(x, function(y) length(E(y))/2)
  })

  Gp.interactome.filt <- Gp.interactome
  for (i in 1:length(Proc.length)){
    Gp.interactome.filt[[i]] <- Gp.interactome.filt[[i]][Proc.length[[i]]>1]
    Proc.length[[i]] <- Proc.length[[i]][Proc.length[[i]]>1]
  }


  annotation.interactome <- list()
  for (i in names(Gp.interactome.filt)){
    annotation.interactome[[i]] <- ann.splitted[[i]][match(names(Gp.interactome.filt[[i]]),
                                                           ann.splitted[[i]]$description),]
  }

  categories.interactome <- lapply(Gp.interactome.filt, function(x){
    prot.categories <- unique(unlist(lapply(x, function(y) names(V(y)))))
    ppi.categories <- induced_subgraph(g.inter.final.hs,prot.categories)
  })

  for (i in names(annotation.interactome)){
    annotation.interactome[[i]]$N.PPI <- Proc.length[[i]]
  }
  return(list(map_id.merged = map_id.merged,
              Gp.interactome = Gp.interactome,
              Gp.interactome.filt = Gp.interactome.filt,
              interactome = interactome.final,
              categories.interactome = categories.interactome,
              annotation.interactome = annotation.interactome))
}
