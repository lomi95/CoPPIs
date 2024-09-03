#' summaryTable
#'
#' @param list_directories list of directories that contains coppis output
#' @param rm_part_categ character vector taht contain eventual part to remove
#' @param rm_part.wb character to remove from worksheet
#' @param name_analysis name to save the summary
#'
#' @importFrom openxlsx createWorkbook
#' @importFrom openxlsx addWorksheet
#' @importFrom openxlsx writeDataTable
#' @importFrom openxlsx saveWorkbook
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#'
#' @return 0 Value if all correct
#' @export

summaryTable <- function(list_directories, rm_part_categ, rm_part.wb,name_analysis){
  wbook <- createWorkbook()
  for (dir.i in list_directories){
    setwd(dir.i)
    xl <- dir()[!(grepl(".RDa",dir()) | grepl("cys",dir()))]
    namesxl <- gsub(paste0(dir.i,"_"), "",dir()[!(grepl(".RDa",dir()) | grepl("cys",dir()))])
    namesxl <- gsub(".xlsx","",namesxl)

    for (i in rm_part_categ){
      namesxl <- gsub(i,"",namesxl)
    }

    names(xl) <- namesxl
    names(namesxl) <- namesxl

    # Aggiunge la colonna categoria
    xl.categ <- lapply(namesxl, function(x){
      xli <- read.xlsx(xl[[x]])
      xli$category <- x
      remaining_cols <- setdiff(colnames(xli), "category")

      # Riorganizza le colonne mettendo 'target_column' in seconda posizione
      xli <- xli[, c(remaining_cols[1], "category", remaining_cols[-1])]
      return(xli)
    })

    df.complete <- Reduce(rbind,xl.categ)
    score.index <- grep("score",colnames(df.complete))
    df.complete.1 <- lapply(1:length(score.index), function(i){

      dd <- df.complete %>% filter(!is.na(df.complete[,score.index[i]]))

      col_to_move <- colnames(dd)[(score.index[i]-1):score.index[i]]
      remaining_cols <- setdiff(colnames(dd), col_to_move)
      dd <- dd[, c(remaining_cols[1:2], col_to_move, remaining_cols[3:length(remaining_cols)])]
      toNULL <- score.index[-1]
      dd[,c(toNULL,toNULL-1)] <- NULL

      addWorksheet(wbook,gsub(rm_part.wb,"",paste0(dir.i,"_",i)))
      writeDataTable(wbook,gsub(rm_part.wb,"",paste0(dir.i,"_",i)), dd)

    })
    setwd("..")
  }
  if (is.null(name_analysis)){
    str_save <- paste0(name_analysis,"_CoPPIs.xlsx")
  } else {
    str_save <- "CoPPIs_summary.xlsx"
  }
  saveWorkbook(wbook, str_save)

}
