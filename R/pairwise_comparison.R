#' Title
#'
#' @param tab_k matrix of k ratios
#' @param g1 index of group 1
#' @param g2 index of group 2
#' @param sd_groups standard deviation of correlations
#' @param tab_ppicor table with PPI and correlation associated
#' @param perc.sign.edges.f table with percentage of significant edges on total edges
#' @param perc.sign filtering threshold, percentage of significant edges on total edges
#' @param significance.coppi significance thresold for CoPPIs pvalue
#' @param correction.coppi correction method for CoPPIs pvalues
#'
#' @importFrom stats pnorm
#' @importFrom stats p.adjust
#'
#' @return A score data.frame
#' @export
#'
#'
#'
pairwise_comparison <- function(tab_k,
                                g1,g2,
                                sd_groups,
                                tab_ppicor,
                                perc.sign.edges.f,
                                perc.sign,
                                significance.coppi,
                                correction.coppi){
  C_cor <- function(Obs,cor1,cor2,tot){
    alpha <- (Obs/tot)/50*log(tot)
    #alpha <- 0.01
    C1 <- (1 - exp(-alpha*(Obs)))/(1 - exp(-alpha*tot))
    return(C1*exp((cor1-cor2)/log(Obs)))
  }

  sd.x <- sd_groups[g1]^2+sd_groups[g2]^2
  if (tab_k[g1,g2] > 1){
    k_x <- 1
    sd.x <- sd.x - (tab_k[g1,g2] - 1)
  } else {
    k_x <- tab_k[g1,g2]
  }


  ratiog1g2 <- cbind(tab_ppicor[,"N.edges"],tab_ppicor[,g1+1]/tab_ppicor[,g2+1])
  if (sum(perc.sign.edges.f[,g2] > perc.sign)){
    ratiog1g2 <- matrix(ratiog1g2[perc.sign.edges.f[,g2] > perc.sign,],
                        ncol = 2)
    if (sum(perc.sign.edges.f[,g2] > perc.sign)==1){
      rownames(ratiog1g2) <- rownames(perc.sign.edges.f)
    } else {
      rownames(ratiog1g2) <- names(which((perc.sign.edges.f[,g2]> perc.sign)==T))
    }

    p.v <- apply(ratiog1g2,1, function(x){
      return(pnorm(x[2],mean = k_x, sqrt(sd.x/x[1])))
    })

    if (is.null(correction.coppi)){
      p.v.adj <- p.v
    } else {
      p.v.adj <- p.adjust(p.v,correction.coppi)
    }
    p.sign <- p.v.adj[p.v.adj < significance.coppi & ratiog1g2[,2] < k_x]



    if (length(p.sign)>0){
      n.ppi <- tab_ppicor[names(p.sign),"N.ppi"]
      n.edg <- tab_ppicor[names(p.sign),"N.edges"]
      z.scr <- (k_x - ratiog1g2[names(p.sign),2])/(-log10(sqrt(sd.x/ratiog1g2[names(p.sign),1])))
      scores1 <- data.frame(description = names(p.sign),
                            N.edges = n.edg,
                            N.PPI   = n.ppi,
                            p.adj   = p.sign,
                            score   = mapply(funscore,n.edg,n.ppi)*sqrt(z.scr)
                            # score  = mapply(C_cor,n.edg,
                            #                    tab_ppicor[names(p.sign),g1+1],tab_ppicor[names(p.sign),g2+1],
                            #                    n.ppi)
      )

    } else {
      scores1 <- data.frame(description = character(),
                            N.edges = numeric(),
                            N.PPI   = numeric(),
                            p.adj   = numeric(),
                            score   = numeric())
    }
  } else {
    scores1 <- data.frame(description = character(),
                          N.edges = numeric(),
                          N.PPI   = numeric(),
                          p.adj   = numeric(),
                          score   = numeric())
  }







  sd.x <- sd_groups[g1]^2+sd_groups[g2]^2
  if (tab_k[g2,g1] > 1){
    k_x <- 1
    if ( (sd.x - (tab_k[g2,g1] - 1))>0 ){
      sd.x <- sd.x - (tab_k[g2,g1] - 1)
    }
  } else {
    k_x <- tab_k[g2,g1]
  }

  ratiog2g1 <- cbind(tab_ppicor[,"N.edges"],tab_ppicor[,g2+1]/tab_ppicor[,g1+1])
  if (sum(perc.sign.edges.f[,g1] > perc.sign)){
    ratiog2g1 <- matrix(ratiog2g1[perc.sign.edges.f[,g1] > perc.sign,],
                        ncol = 2)

    if (sum(perc.sign.edges.f[,g1] > perc.sign)==1){
      rownames(ratiog2g1) <- rownames(perc.sign.edges.f)[perc.sign.edges.f[,g1] > perc.sign]
    } else {
      rownames(ratiog2g1) <- names(which((perc.sign.edges.f[,g1]> perc.sign)==T))
    }

    p.v <- apply(ratiog2g1,1, function(x){
      return(pnorm(x[2],mean = k_x, sqrt(sd.x/x[1])))
    })

    if (is.null(correction.coppi)){
      p.v.adj <- p.v
    } else {
      p.v.adj <- p.adjust(p.v,correction.coppi)
    }
    p.sign <- p.v.adj[p.v.adj < significance.coppi & ratiog2g1[,2] < k_x]

    if (length(p.sign)>0){
      n.ppi <- tab_ppicor[names(p.sign),"N.ppi"]
      n.edg <- tab_ppicor[names(p.sign),"N.edges"]
      z.scr <- (k_x - ratiog2g1[names(p.sign),2])/(-log10(sqrt(sd.x/ratiog2g1[names(p.sign),1])))
      scores2 <- data.frame(description = names(p.sign),
                            N.edges = n.edg,
                            N.PPI   = n.ppi,
                            p.adj   = p.sign,
                            score   = -mapply(funscore,n.edg,n.ppi)*sqrt(z.scr)
                            # score  = -mapply(C_cor,n.edg,
                            #                  tab_ppicor[names(p.sign),g2+1], tab_ppicor[names(p.sign),g1+1],
                            #                  n.ppi)
      )


    } else {
      scores2 <- data.frame(description = character(),
                            N.edges = numeric(),
                            N.PPI   = numeric(),
                            p.adj   = numeric(),
                            score   = numeric())
    }
  } else {
    scores2 <- data.frame(description = character(),
                          N.edges = numeric(),
                          N.PPI   = numeric(),
                          p.adj   = numeric(),
                          score   = numeric())
  }


  scores <- rbind(scores1,scores2)
  return(scores)
}
