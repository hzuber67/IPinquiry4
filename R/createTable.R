#' Interactive table
#'
#' \code{createTable} create an html interactive table using DT package
#'
#' @param tres Table containing an header and one name for each lines
#' @return Return interactive table in html format
#' @details
#' Interactive table containing a search option and allowing to select and export slected data using csv and xls.
#' @export
#'
createTable <- function(tres)
{
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("Package \"DT\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  DT::datatable(tres,
            extensions = c('FixedHeader','ColReorder', 'Buttons'),
            options = list(pageLength = 50, fixedHeader = TRUE, colReorder = TRUE,
                           dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel'
                           )))
}
