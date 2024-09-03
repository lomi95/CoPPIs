#' Title
#'
#' @param file_path directory path in which the CoPPIs output excels are
#' @param filter_similarity filtering threshold, minimum similarity score between terms
#'
#' @importFrom openxlsx getSheetNames
#' @importFrom openxlsx read.xlsx
#' @importFrom RCy3 createNetworkFromDataFrames
#' @importFrom RCy3 loadTableData
#' @return message of loading states
#' @export
#'
Load2Cytoscape <- function(file_path, filter_similarity = 0.6){
  all_sheets <- getSheetNames(file_path)
  names(all_sheets) <- all_sheets
  for (i in all_sheets[-1]){
    message("__________________________________________\n",file_path," - ",i)
    xl.i <- read.xlsx(file_path,sheet = i)
    nodes <- unique(data.frame(id = xl.i[,1],
                               scores = as.numeric(xl.i[,"score1"])),
                    stringsAsFactors=FALSE)
    edges <- data.frame(source=xl.i[,1],
                        target=xl.i[,2],
                        weight=as.numeric(xl.i[,3]),
                        stringsAsFactors=FALSE)
    if (nrow(nodes)){
      createNetworkFromDataFrames(nodes,edges[edges$weight > filter_similarity,],
                                  title= i, collection= gsub(".xlsx","",file_path))
      loadTableData(read.xlsx(file_path,all_sheets[1]), data.key.column = "description")
    }
  }
}
