#' Function to calculate the scale factor as done by DEseq package (Median-to-ratio method)
#'
#' Function to calculate the scale factor using the median to ratio method as used in DEseq2
#' @param data Objet de type IPObj, provenant directement de la fonction \code{loadData}.
#' @return Return the scale factor calculated using the median to ratio method as used in DEseq2
#' @importFrom stats median
#
#
median_of_ratios <- function(data)
{   gmean <- function(x, na.rm = FALSE){
  if(na.rm) x <- x[!is.na(x)]
  n <- length(x)
  prod(x)^(1/n)
}
div <- apply((data+1)/ apply(data + 1,1,gmean),2, median)
div
}
