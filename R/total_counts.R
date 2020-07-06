#' Fonction to count the total number of counts for each IP
#' Helene Zuber
#' \code{total_counts} This function allow for calculating the number of spectra for each IP
#' This could then be used as scale factor for statistical analysis
#'
#' @param data Objet de type IPObj, that comes from the function \code{loadData}.
#' @return The function return the total number of counts by sample
#'
#'
total_counts <- function(data)
{   div<- apply(data,2,sum)
div
}
