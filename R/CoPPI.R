#' Title
#'
#' @param annotations special annotation object
#' @param correlations special correlation object
#' @param perc.ned.nppi filtering threshold, minimum percentage of significant number of edges divided by number of PPIs
#' @param perc.sign filtering threshold, minimum percentage of significant number of edges on total number of edges
#' @param significance.coppi significance threshold for CoPPIs terms pvalues
#' @param correction.coppi correction method for CoPPIs terms
#' @param min_edges filtering threshold, minimum edges of a biological term
#' @param max_edges filtering threshold, maximum edges of a biological term
#' @param names_of_groups names of the groups
#' @param p.filter filtering threshold, maximum CoPPIs pvalue
#' @param score.filter filtering threshold, minimum CoPPIs score
#' @param interactome interactome as data.frame
#'
#'
#' @importFrom igraph delete_vertices
#' @importFrom igraph degree
#' @importFrom igraph E
#' @importFrom igraph make_empty_graph
#' @importFrom igraph subgraph
#' @importFrom igraph V
#' @importFrom igraph vertex.attributes
#' @importFrom igraph vertex.attributes<-
#' @importFrom matrixStats rowMaxs
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#'
#' @importFrom dplyr filter
#' @importFrom dplyr %>%
#'
#'
#' @return List with information about CoPPIs result
#' @export
#'
#'

CoPPI <- function(annotations,
                  correlations,
                  perc.ned.nppi,
                  perc.sign,
                  significance.coppi,
                  correction.coppi,
                  min_edges,
                  max_edges,
                  names_of_groups,
                  p.filter,
                  score.filter,
                  interactome){
  st.time <- Sys.time()
  Other.PPI <- correlations$Other.PPI
  sd_groups <- correlations$Sd_groups
  proteins <- unique(unlist(annotations$genes_found))
  pathways <- annotations$description
  tab_pathways_protein <- matrix(ncol=length(proteins),nrow = length(pathways),data = 0,
                                 dimnames = list(pathways,proteins))

  for (i in 1:nrow(annotations)){
    tab_pathways_protein[i,match(unlist(annotations[i,]$genes_found),proteins)] <- 1
  }
  #creating path-path network
  tab_protein_pathways <- t(tab_pathways_protein)

  flattCorr <- correlations$flattCorr


  gPath.All <- lapply(correlations$gCor.categ.all, function(x){
    apply(tab_protein_pathways,2,function(y){
      vids <- intersect(names(which(y==1)),names(V(x)))
      subG.0 <- subgraph(x,vids)
      subG <- delete_vertices(subG.0, which(degree(subG.0)==0))
      if (length(E(subG)) >= min_edges && length(E(subG)) <= max_edges){
        return(subG)
      } else {
        empty <- make_empty_graph(directed = F)
        return(empty)
      }
    })
  })

  gPath.Sign <- lapply(correlations$gCor.categ.sign, function(x){
    apply(tab_protein_pathways,2,function(y){
      vids <- intersect(names(which(y==1)),names(V(x)))
      subG.0 <- subgraph(x,vids)
      subG <- delete_vertices(subG.0, which(degree(subG.0)==0))
      if (length(E(subG)) >= min_edges && length(E(subG)) <= max_edges){
        return(subG)
      } else {
        empty <- make_empty_graph(directed = F)
        return(empty)
      }
    })
  })

  # TAB K #
  ######
  N.g <- length(names_of_groups)
  # calcoliamo il rapporto tra la media dei pesi dei gruppi (k)
  tab_k <- matrix(nrow=length(flattCorr),ncol=length(flattCorr))
  for (i in 1:nrow(tab_k)){
    for (j in 1:ncol(tab_k)){
      tab_k[i,j] <- mean(c(abs(flattCorr[[i]]$All_correlation$weights),rep(1,Other.PPI)))/
        mean(c(abs(flattCorr[[j]]$All_correlation$weights),rep(1,Other.PPI)))
    }
  }
  colnames(tab_k) <- names(flattCorr)
  rownames(tab_k) <- names(flattCorr)


  if (length(gPath.All[[1]])){

    tab_meancor <- matrix(nrow=length(gPath.All[[1]]),
                          ncol=length(flattCorr)+2,
                          dimnames = list(names(gPath.All[[1]]),
                                          c("N.edges",
                                            paste("cor.",names(flattCorr)),
                                            "N.ppi")))
    for (i in 1:length(flattCorr)){
      e.i <- unlist(mean_abs_cor(gPath.All,i))
      tab_meancor[names(e.i),1]   <- unlist(sapply(gPath.All[[i]],function(x){
        length(E(x))
      }))
      tab_meancor[names(e.i),i+1] <- e.i
    }
    tab_meancor[,ncol(tab_meancor)] <- annotations$N.PPI
    tab_meancor <- tab_meancor[order(tab_meancor[,1]),]
    tab_meancor <- tab_meancor[!is.na(tab_meancor[,1]),]

    ## si calcola la media considerando le correlazioni == 1 per tutte le ppi non osservate
    tab_ppicor0 <- t(apply(tab_meancor,1,function(x){
      x[2:(length(x)-1)] <- (x[2:(length(x)-1)]*x[1]+(x[length(x)]-x[1]))/x[length(x)]
      return(x)
    }))

    tab_ppicor1 <- tab_ppicor0[(tab_ppicor0[,1]/tab_ppicor0[,ncol(tab_ppicor0)])>perc.ned.nppi,]

    perc.sign.edges.f <- sapply(gPath.Sign, function(x){
      sapply(x[rownames(tab_ppicor1)], function(y) length(E(y)))
    })/tab_ppicor1[,1]

    if (sum(rowMaxs(perc.sign.edges.f)>perc.sign)){
      tab_ppicor        <- matrix(tab_ppicor1[rowMaxs(perc.sign.edges.f)>perc.sign,],
                                  ncol = ncol(tab_ppicor1),
                                  dimnames = list(names(which(rowMaxs(perc.sign.edges.f)>perc.sign)),
                                                  colnames(tab_ppicor1)))
      perc.sign.edges.f <- matrix(perc.sign.edges.f[rowMaxs(perc.sign.edges.f)>perc.sign,],
                                  ncol = ncol(perc.sign.edges.f),
                                  dimnames = list(names(which(rowMaxs(perc.sign.edges.f)>perc.sign)),
                                                  colnames(perc.sign.edges.f)))
      # funzione che prende i significativi e assegna lo score
      sig_path <- list()
      n <-1
      for (i in 1:(length(flattCorr)-1)){
        for (j in (i+1):length(flattCorr)){

          sig_path[[n]] <- pairwise_comparison(tab_k = tab_k,
                                               g1 = i,g2 = j,
                                               sd_groups = sd_groups,
                                               tab_ppicor = tab_ppicor,
                                               perc.sign.edges.f = perc.sign.edges.f,
                                               perc.sign = perc.sign,
                                               significance.coppi = significance.coppi,
                                               correction.coppi = correction.coppi)

          sig_path[[n]] <- sig_path[[n]] %>%
            filter(p.adj < p.filter &
                     abs(score) > score.filter &
                     N.edges/N.PPI > perc.sign)
          tryCatch({
            names(sig_path)[n] <- paste(names(flattCorr)[j]," vs ",names(flattCorr)[i])
            n <- n+1
          },error = function(err){
            if (length(table(annotations$category))==1){
              message("No significant terms were found in ",
                      paste(names(flattCorr)[j]," vs ",names(flattCorr)[i]),
                      " in ",annotations$category[1])
            } else {
              message("No significant terms were found in ",
                      paste(names(flattCorr)[j]," vs ",names(flattCorr)[i]))
            }
          })
        }
      }

      ## in tabella: annotations con risultati
      all_processes <- sort(unique(unlist(sapply(sig_path, function(x) x$description))))
      if (length(all_processes)){
        final_table.0 <- annotations[match(all_processes,annotations$description),]
        if (nrow(tab_ppicor)>1){
          final_table.1 <- cbind.data.frame(final_table.0,
                                            tab_ppicor[final_table.0$description,])
        } else {
          final_table.1 <- cbind.data.frame(final_table.0,tab_ppicor)
        }
        final_table.1$N.ppi <- NULL

        sign.edges <- sapply(gPath.Sign, function(x){
          sapply(x[final_table.1$description], function(y) length(E(y)))
        })
        if (is.null(dim(sign.edges))){
          names.s <- names(sign.edges)
          sign.edges <- matrix(sign.edges, nrow = 1)
          rownames(sign.edges) <- unique(sapply(str_split(names.s, "\\."), function(x) x[2]))
          colnames(sign.edges) <- sapply(str_split(names.s, "\\."), function(x) x[1])
        }
        perc.edges <- sign.edges/final_table.1$N.edges

        colnames(sign.edges) <- paste("Sign_edges.",colnames(sign.edges))
        colnames(perc.edges) <- paste("Perc_edges.",colnames(perc.edges))


        Nodes.PPI <- sapply(gPath.All[[1]][final_table.1$description], function(x){
          length(V(x))
        })
        sign.nodes <- sapply(gPath.Sign, function(x){
          sapply(x[final_table.1$description], function(y) length(V(y)))
        })
        if (is.null(dim(sign.nodes))){
          names.n <- names(sign.nodes)
          sign.nodes <- matrix(sign.nodes, nrow = 1)
          rownames(sign.nodes) <- unique(sapply(str_split(names.s, "\\."), function(x) x[2]))
          colnames(sign.nodes) <- sapply(str_split(names.s, "\\."), function(x) x[1])
        }
        perc.nodes <- sign.nodes/Nodes.PPI
        colnames(sign.nodes) <- paste("Sign_nodes.",colnames(sign.nodes))
        colnames(perc.nodes) <- paste("Perc_nodes.",colnames(perc.nodes))

        final_table.2 <- cbind.data.frame(final_table.1,
                                          Nodes.PPI,
                                          sign.nodes,
                                          sign.edges,
                                          perc.nodes,
                                          perc.edges)


        scorespvalue <- lapply(sig_path,function(x){
          x[,c("p.adj","score")]
        })
        final.table <- final_table.2
        for (i in names(scorespvalue)){
          colnames(scorespvalue[[i]]) <- paste(colnames(scorespvalue[[i]])," ",i)
          scorespvalue[[i]]$description <- rownames(scorespvalue[[i]])
          final.table <- merge(final.table,scorespvalue[[i]],
                               by = "description", all=T)
        }
        final.table <- final.table[!is.na(final.table$term),]



        graph.similarity <- lapply(sig_path, function(y){
          if (sum(is.na(y))){
            graph_terms <- make_empty_graph(directed = F)
            el_gpp <- data.frame(Node_path1 = NA,
                                 Node_path2 = NA,
                                 Similarity = NA,
                                 p.adj1 = NA,
                                 score1 = NA,
                                 p.adj2 = NA,
                                 score2 = NA)
            return(list(graph_term = graph_terms,
                        el_graph.term = el_gpp))
          } else {
            W <- tab_pathways_protein[y$description,]
            if (is.null(dim(W))){
              el_gpp <- data.frame(Node_path1 = y$description,
                                   Node_path2 = y$description,
                                   Similarity = 1,
                                   p.adj1 = y$p.adj,
                                   score1 = y$score,
                                   p.adj2 = y$p.adj,
                                   score2 = y$score)
              graph_terms <- graph_from_edgelist(as.matrix(el_gpp[,1:2]))
              vertex.attributes(graph_terms)$scores <- y$score
              return(list(graph_term = graph_terms,
                          el_graph.term = el_gpp))
            } else {
              W <- W[,colSums(W)!=0]
              return(graph_terms(W,y))
            }
          }
        })

      } else {
        final.table <- NULL
        graph.similarity <- NULL
      }
    } else {
      sig_path      <- NULL
      final.table   <- NULL
      graph.similarity <- NULL
    }



  } else {
    sig_path      <- NULL
    final.table   <- NULL
    graph.similarity <- NULL
  }

  if (length(table(annotations$category))==1){
    time_employed <- format_time_difference(difftime(Sys.time(),st.time,units = "s"))
    message(paste(annotations$category[1],"finished in", time_employed,sep = " "))
  }



  return(list(gPath.All  = gPath.All,
              gPath.Sign = gPath.Sign,
              sig_path   = sig_path,
              tab_pathways_protein = tab_pathways_protein,
              graph.similarity = graph.similarity,
              table_summary = final.table))

}
