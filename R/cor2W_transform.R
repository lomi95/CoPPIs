#' Title
#'
#' @param flattCorr correlation martix flatten as data.frame
#' @param signifCorr significance threshold
#'
#' @return correlation matrix flatten as data.frame with weights column
#' @export
#'
#'
cor2W_transform <- function(flattCorr, signifCorr){

  ### normalizzazione per togliere gli 1 dal pvalue, per rendere possibile il logaritmo
  # L'obiettivo è trasformare l'intervallo del pvalue da (Significativià ; 1 ] a (0 ; 1)

  # Il nuovo massimo del pvalue sarà la media tra 1 e il secondo massimo
  if (max(flattCorr$p.adj)==1){
    new_max <- mean(c(1,max(flattCorr$p.adj[-which(flattCorr$p.adj==1)])))
  } else {
    new_max <- max(flattCorr$p.adj)
  }

  # calcolo parametri normalizzazione: (pvalue + epsilon ) * k = p_norm1
  # ( signifCorr + eps ) * k = signifCorr
  # ( max(pvalue)  + eps ) * k = new_max

  # risolvendo questo sistema si ottiene
  eps <- ( new_max - max(flattCorr$p.adj) ) * signifCorr / ( signifCorr - new_max )
  k   <- ( signifCorr - new_max ) / ( signifCorr - max(flattCorr$p.adj) )

  # applichiamo la normalizzazione solo ai pvalue non significativi
  ind_notsig <- which(flattCorr$p.adj > signifCorr)
  norm1 <- (flattCorr$p.adj[ind_notsig] + eps) * k

  # adesso applichiamo -log10 per "invertire" i pvalue bassi in pesi alti e viceversa
  lg    <- -log10(norm1)

  ## normalizziamo lg tra 0 e 1, calcolo parametri ( lg + eps_log) * k_log = p_norm2
  # ( MAX_lg  + eps_log ) * k_log = 1
  # ( min(lg) + eps_log ) * k_log = min(lg)
  MAX_lg <- -log10(signifCorr)

  # risolvendo il sistema si ottiene
  eps_log <- ( min(lg) - MAX_lg * min(lg) ) / ( min(lg) - 1 )
  k_log   <- ( min(lg) - 1 ) / ( min(lg) - MAX_lg)
  norm2   <- (lg + eps_log) * k_log

  #moltiplichiamo la correlazione per il pnorm3
  flattCorr$weights <- flattCorr$cor
  flattCorr$weights[ind_notsig] <- flattCorr$cor[ind_notsig]*norm2^3

  return(flattCorr)
}
