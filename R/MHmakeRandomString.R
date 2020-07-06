#' Create sub dataset
#' Timothee Vincent
#' Generated 12 letter random string
#'
#' \code{MHmakeRandomString} Annexe function called by the function \code{moreThantwoCond}. It is used to create temprory file.
#' This file will be erased at the end
#' @param n (default 1). Number of generated randomstring
#' @param length (default 12), Size of the random string
#' @return Return a random string

MHmakeRandomString <- function(n=1, length=12)
{
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    length, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}
