#' Create Multidimensional scaling (MDS) plots
#'
#' Function to plot data in a 2D plan such that the between-object distances are preserved as well as possible.
#' \code{MDSplot} is  useful to visualise the level of similarity between experiments
#'
#' @param IPObj object that comes directly from the function \code{load_IP_Data}
#' @param norm (default nothing), if \code{norm}="total", all counts are divided by the total number of count,
#' if \code{norm}="DEseq", counts are normalized according to the median-to-ratios methods as used in the DEseq2 R package
#' @return Multidimensional scaling (MDS) plot based on euclidean distances
#' @importFrom stats dist cmdscale
#' @importFrom graphics plot text abline
#' @export
MDSplot <- function(IPObj, norm = c("nothing","total","DEseq"))
{
  norm <- match.arg(norm)
  data <- (IPObj$spectral_counts)
  if (norm=="nothing")
  {
    trans <- t(data)
    d <- dist(trans) # to calculate euclidean distances between the rows
    fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
    normInfo <- "No normalisation"
  }
  else if (norm=="total")
  {
    norm <- t(t(data)/total_counts(data)) # all counts are divided by the total number of count per sample
    trans <- t(norm)
    d <- dist(trans) # euclidean distances between the rows
    fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dimensions
    normInfo <- "Normalisation based on the total count number"
  }
  else if (norm=="DEseq")
  {
    norm <- t(t(data)/median_of_ratios(data)) * mean(median_of_ratios(data))
    trans <- t(norm)
    d <- dist(trans) # euclidean distances between the rows
    fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
    normInfo <- "Median-to-ratio normalisation"
  }

  pcEig <- fit$eig * 100 / sum(fit$eig)
  PC1 <- round(pcEig[1], 1)
  PC2 <- round(pcEig[2], 1)
  x <- fit$points[,1]
  y <- fit$points[,2]
  x_max <- max(x) + (abs(max(x)) * 0.3)
  y_max <- max(y) + (abs(max(y)) * 0.3)
  x_min <- min(x) - (abs(min(x)) * 0.3)
  y_min <- min(y) - (abs(min(y)) * 0.3)
  plot(x, y, xlab=sprintf("Coordinate 1 (%.2f%%)", PC1), ylab=sprintf("Coordinate 2 (%.2f%%)", PC2), main="Metric MDS", sub=sprintf("%s", normInfo),
       type="p",pch = 20, asp = 1, cex.main=1.2, cex.lab=1.0, cex.sub=1.0, font.main = 2, font.sub = 3,col=c(IPObj$treat), xlim=c(x_min,x_max), ylim=c(y_min,y_max))
  text(x, y, labels = substr(row.names(trans), 1, 20), cex=0.7,col=c(IPObj$treat),font=1, offset = 0.5, adj = c(0.5,2))
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
}
