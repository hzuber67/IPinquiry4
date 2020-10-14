#' Add supplemental information column(s) to the result table
#'
#' \code{add_suppl_information} generates a new table with new information in supplemental columns

#' @param table table with statistical results with protein or gene ID as row.names
#' @param info .txt file with at least two columns and with an header. first column correspond to gene of protein ID
#' other columns correspond to interest information E
#' @return Combine \code{table} and \code{info} table
#' @export
#'
add_suppl_information <- function(table, info)
{
  supp <- read.table(info, header=T, sep= "\t", stringsAsFactors=FALSE, quote="", fill=FALSE) ##table with list of AGI for RH protein
  row.names(supp) <- supp[,1]
  supp <- supp[,2:ncol(supp),drop=F]
  mergeObj <- merge(table, supp, by="row.names", all.x=TRUE)
  row.names(mergeObj) <- mergeObj[,1]
  mergeObj <- mergeObj[,2:ncol(mergeObj)]
  mergeObj
}
