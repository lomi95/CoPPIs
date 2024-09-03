#' graph terms
#'
#' @param W Weights
#' @param S Scores
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph vertex.attributes
#' @importFrom igraph edge.attributes
#' @importFrom igraph as_edgelist
#' @importFrom igraph vertex.attributes<-
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats colMaxs
#'
#' @return List information about terms similarity network
#' @export
#'
graph_terms <- function(W,S){

  if (is.null(dim(W))){
    W1 <- t(W) %*% W
    colnames(W1) <- S$description
    rownames(W1) <- S$description

  } else {
    W1 <- W %*% t(W)
  }
  W2 <- (W1/rowMaxs(W1) + t(W1/colMaxs(W1)))/2



  graph_terms <- graph_from_adjacency_matrix(W2,mode="undirected",
                                             weighted = T, diag = T)
  vertex.attributes(graph_terms)$scores <- as.numeric(S$score)
  el <- as_edgelist(graph_terms)
  el_gpp <- data.frame(Node_path1 = el[,1],
                       Node_path2 = el[,2],
                       Similarity = edge.attributes(graph_terms)$weight,
                       p.adj1 = S[el[,1],"p.adj"],
                       score1 = S[el[,1],"score"],
                       p.adj2 = S[el[,1],"p.adj"],
                       p.adj1 = S[el[,1],"score"])

  return(list(graph_term = graph_terms,
              el_graph.term = el_gpp))
}
