#' Title
#'
#' @param list.groups list of dataset
#' @param names_of_groups names of the groups
#' @param categories.interactome interactome for each cateogry
#' @param gCOR.groups.l if a precedent analysis was run, an object gCOR.groups.l may be used to speed up the computation
#' @param min_n.corr minimum number of observation pairs for computing the correlation
#' @param signifCorr significance threshold for correlation pvalue
#' @param correctionCorr pvalue adjusting method
#' @param corr.test type of correlation test (Spearman or Pearson)ù
#' @param gInter graphs of interactions
#' @param compute_weights boolean, the weight transformation to the non significant
#'     correlation need to be computed?
#'
#' @importFrom igraph edge.attributes<-
#' @importFrom stats p.adjust
#' @importFrom stats sd
#'
#'
#'
#' @return a list containing information about correlation of conditions
#' @export
#'
Correlations_CoPPIs <- function(list.groups,
                                names_of_groups,
                                categories.interactome,
                                gCOR.groups.l,
                                min_n.corr,
                                signifCorr,
                                correctionCorr,
                                corr.test,
                                gInter,
                                compute_weights
){
  start.time.correlation <- Sys.time()
  CorrNULL <- matrix(1, nrow = ncol(list.groups[[1]]), ncol = ncol(list.groups[[1]]),
                     dimnames = list(colnames(list.groups[[1]]),colnames(list.groups[[1]])))
  flattNULL <- flattenCorrMatrix(CorrNULL,CorrNULL,CorrNULL)
  flattNULL$cor_features <- paste(flattNULL$row,"and",flattNULL$column)

  # intersection with interactome of categories
  gNULL.categories   <- lapply(categories.interactome, function(y){

    gNULL <- graph_from_edgelist(as.matrix(flattNULL[,1:2]),directed = F)
    edge.attributes(gNULL)$features <- flattNULL$cor_features

    gNULL.inters <- intersection(gNULL, y, keep.all.vertices = F)
    ind.ppi <- match(sort(edge.attributes(gNULL.inters)$features),flattNULL$cor_features)
    gNULL.el <- flattNULL[ind.ppi,]

    return(list(el = gNULL.el,
                length.interactome = length(E(y))))
  })
  names(gNULL.categories) <- str_to_title(names(gNULL.categories))

  if (is.null(gCOR.groups.l)){
    gCOR.groups.l <- lapply(list.groups,function(x){
      Corr <- rcorr(x,type = corr.test)
      flattCorr <- flattenCorrMatrix(Corr$r,Corr$P,Corr$n)
      flattCorr$cor_features <- paste(flattCorr$row,"and",flattCorr$column)
      flattCorr$cor[is.na(flattCorr$cor)] <- 0
      flattCorr$p[is.na(flattCorr$p)] <- 1

      if (!is.null(min_n.corr)){
        flattCorr$p[flattCorr$n < min_n.corr,] <- 1
        compute_weights <- T
      }

      if (is.null(correctionCorr)){
        flattCorr$p.adj <- flattCorr$p
      } else {
        flattCorr$p.adj <- p.adjust(flattCorr$p, correctionCorr)
      }
      ind_sig <- flattCorr$p.adj <= signifCorr

      if (compute_weights){
        flattCorr <- cor2W_transform(flattCorr,signifCorr)

      } else {
        flattCorr$weights <- flattCorr$cor
        flattCorr$weights[!ind_sig] <- 0
      }

      Sign <- flattCorr[ind_sig,] #storing significative ones
      Sign <- Sign[order(Sign$cor_features),]

      netw.weighted <- flattCorr[,c("row","column","weights")]
      gcor <- graph_from_edgelist(as.matrix(netw.weighted[,1:2]),directed = F)
      edge.attributes(gcor)$weights <- netw.weighted[,3]

      g.f <- intersection(gcor,gInter, keep.all.vertices = F) ## grafo correlazioni e interazioni
      return(list(gCor.groups   = flattCorr,
                  gInter.groups = g.f))

    })
  }
  gCOR.groups   <- lapply(gCOR.groups.l, function(x) x$gCor.groups)

  gCOR.categories <- lapply(gNULL.categories, function(x){
    ## aggiungere trycatch con mesasggio di funzione gcor.groups standard
    flattCateg <- lapply(gCOR.groups, function(y){
      gcor.el <- y[match(x$el$cor_features,y$cor_features),]
      if (is.null(correctionCorr)){
        gcor.el$p.adj <- gcor.el$p
      } else {
        gcor.el$p.adj <- p.adjust(gcor.el$p, correctionCorr)
      }
      ind_sig <- gcor.el$p.adj <= signifCorr

      if (compute_weights){
        cor2W_transform(gcor.el, signifCorr)

      } else {
        gcor.el$weights <- gcor.el$cor
        gcor.el$weights[!ind_sig] <- 0
      }

      Sign <- gcor.el[ind_sig,] #storing significative ones
      Sign <- Sign[order(Sign$cor_features),]

      return(list(All_correlation  = gcor.el,
                  Sign_correlation = Sign))
    })
    gCor.categ.all <- lapply(flattCateg, function(y){
      g.Cor <- graph_from_edgelist(as.matrix(y$All_correlation[,1:2]), directed = F)
      edge.attributes(g.Cor)$weight <- y$All_correlation$weights
      return(g.Cor)
    })
    gCor.categ.sign <- lapply(flattCateg, function(y){
      g.Cor <- graph_from_edgelist(as.matrix(y$Sign_correlation[,1:2]), directed = F)
      edge.attributes(g.Cor)$weight <- y$Sign_correlation$weights
      return(g.Cor)
    })

    Other.PPI <- x$length.interactome - nrow(x$el)
    sd_groups <- sapply(flattCateg, function(x){
      sd(c(x$All_correlation$weights,rep(1,Other.PPI)))/2
    })
    return(list(flattCorr = flattCateg,
                gCor.categ.all = gCor.categ.all,
                gCor.categ.sign = gCor.categ.sign,
                Sd_groups = sd_groups,
                Other.PPI = Other.PPI))
  })
  time_employed <- format_time_difference(difftime(Sys.time(),
                                                   start.time.correlation,units = "s"))
  message("Correlations computed in ", time_employed,sep = " ")

  return(list(gCOR.groups.l = gCOR.groups.l,
              gCOR.categories = gCOR.categories))
}
