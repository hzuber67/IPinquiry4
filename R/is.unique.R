#' Create sub dataset
#' Timothee Vincent
#' Check that a vector is unique
#'
#' \code{is.unique} use vector as argument and check that all vector elements are unique.
#' If true, return \code{TRUE}, otherwise return \code{FALSE}
#' @param vect a vector
#' @return return \code{TRUE} if the vector contains only unique elements, return \code{FALSE} otherwise
is.unique <- function(vect)
{
  if (length(vect) == length(unique(vect))){TRUE}
  else {FALSE}
}

