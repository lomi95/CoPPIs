#' Title
#'
#' @param SUID_list vector of network SUIDs, default getNetworkList(getSUIDs = T)
#' @param style.name name of the style, default "score"
#' @param palette.nodes palette of the nodes, default c('#4EB3D3','#D4FFD0','#FFE0FF','#FC0000')
#'
#' @importFrom RCy3 getCollectionSuid
#' @importFrom RCy3 getVisualPropertyDefault
#' @importFrom RCy3 mapVisualProperty
#' @importFrom RCy3 createVisualStyle
#' @importFrom RCy3 lockNodeDimensions
#' @importFrom RCy3 getNetworkName
#' @importFrom RCy3 getCollectionName
#' @importFrom RCy3 getCollectionSuid
#' @importFrom RCy3 setVisualStyle
#' @importFrom RCy3 getNetworkList
#' @importFrom RCy3 getVisualPropertyNames
#'
#' @export
#'
#'
#'
setCytoStyle <- function(SUID_list = getNetworkList(getSUIDs = T),
                         style.name = "score",
                         palette.nodes = c('#4EB3D3','#D4FFD0','#FFE0FF','#FC0000')){
  all_Suid.collection <- unique(sapply(SUID_list,getCollectionSuid))
  default_param <- sapply(getVisualPropertyNames(),getVisualPropertyDefault,"default")
  nodeLabels <- mapVisualProperty('node label','description','p') #description was id
  nodeFills <- mapVisualProperty(visual.prop = 'node fill color',table.column = 'scores',mapping.type = 'c',
                                 table.column.values = as.numeric(c(-8,-.Machine$double.eps,
                                                                    .Machine$double.eps,8)),
                                 visual.prop.values = palette.nodes)
  createVisualStyle(style.name, default_param, list(nodeLabels,nodeFills))
  lockNodeDimensions(FALSE, style.name)
  lapply(SUID_list,function(x){
    message(getNetworkName(x), " - ", getCollectionName(getCollectionSuid(x)))
    setVisualStyle(style.name = style.name,network = x)
  })
}
