#' merge_annotations
#'
#' @param df1 first annotation data.frame
#' @param df2 second annotation data.frame
#'
#' @importFrom stats na.omit
#'
#' @return merged annotation data.frame
#' @export
#'
merge_annotations <- function(df1,df2){
  df.m <- merge(df1,df2, by = "term",all = T)

  col.original <- c("category","description")
  col.genes    <- c("inputGenes","preferredNames")

  new_col <- list()
  for (i in col.original){
    col.i <- df.m[,grep(i,colnames(df.m))]
    new_col[[i]] <- apply(col.i,1,function(x){
      x[which(!is.na(x))[1]]
    })
  }
  union.not.na <- function(x,y){
    as.vector(na.omit(union(x,y)))
  }
  for (i in col.genes){
    col.i <- df.m[,grep(i,colnames(df.m))]
    new_col[[i]] <- mapply(union.not.na, col.i[,1], col.i[,2])
  }
  new_col.df <- data.frame(term = df.m$term,
                           category = new_col$category,
                           description = new_col$description)
  new_col.df$inputGenes <- new_col$inputGenes
  new_col.df$preferredNames <- new_col$preferredNames
  new_col.df$number_of_genes <- sapply(new_col$preferredNames,length)
  return(new_col.df)
}
