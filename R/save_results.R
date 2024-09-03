#' Title
#'
#' @param category category object of CoPPIs output
#' @param name_analysis name of the analysis to be saved
#'
#' @importFrom openxlsx addWorksheet
#' @importFrom openxlsx writeDataTable
#' @importFrom openxlsx createWorkbook
#' @importFrom openxlsx saveWorkbook
#'
#'
#' @return message of saving state
#' @export
#'
save_results <- function(category, name_analysis){
  if (!is.null(category$table_summary)){

    name_categ <- category$table_summary$category[1] %>%
      substr(1,15)
    if (!is.null(name_categ)){
      category$table_summary$category <- NULL
      category$table_summary$inputGenes <- NULL
      summary.terms <- createWorkbook()
      # nell'enrichment deve esserci
      # uno sheet con la tabella con tutti i pathways risultati significativi
      name.sheet <- "node_table.terms"

      addWorksheet(summary.terms, name.sheet)
      writeDataTable(summary.terms, name.sheet,
                     category$table_summary)

      for (i in 1:length(category$graph.similarity)){
        if (nrow(category$graph.similarity[[i]]$el_graph.term)){
          name.sheet <- names(category$graph.similarity)[[i]]
          addWorksheet(summary.terms, name.sheet)
          writeDataTable(summary.terms, name.sheet,
                         category$graph.similarity[[i]]$el_graph.term)
        }
      }

      saveWorkbook(summary.terms,
                   paste0(name_analysis,"_",
                         name_categ,".xlsx"),
                   overwrite = TRUE)
    }
  }

}
