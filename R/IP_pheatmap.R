#' Heatmap using the package pheatmap
#'
#' \code{IP_pheatmap} uses the pheatmap package to plot heatmaps showing spectral counts across samples for selected genes/proteins.
#' @param IPObj Object generated from the function \code{loadData}
#' @param GeneList List of selected genes/proteins used to plot for the heatmap.
#' @param norm If \code{"nothing"}, no normalization, if  \code{"total"}, counts are divided by the total number of counts,
#' if \code{"DEseq"}, counts are normalized using the median of ratios (DEseq method)
#' @param cluster_rows Determine if rows should be clustered, it can be \code{TRUE} or \code{FALSE}
#' @param cluster_cols Determine if columns should be clustered, it can be \code{TRUE} or \code{FALSE}
#' @param annotation_col Data frame that specifies the annotations shown on left side of the heatmap.
#' @param annotation_row Data frame that specifies the annotations shown on top of the heatmap.
#' @param fontsize_row To change font size of the row label. By defaut, the font size is 6
#' @param fontsize_col To change font size of the col label. By defaut, the font size is 6
#' @param title To change graph title
#' @details
#' This function extract count data for selected genes \code{GeneList} and create an heatmap for all IP in the \code{IPObj}
#' @return Return heatmap for selected genes/proteins
#' @export
#'
IP_pheatmap <- function(IPObj, GeneList, norm = c("nothing", "total", "DEseq"), cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=NULL, annotation_row=NULL, fontsize_row=6, fontsize_col=6, title=NULL)
{
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package \"pheatmap\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  norm <- match.arg(norm)
  data <- (IPObj$spectral_counts)
  if (norm=="nothing")
  { MyCounts <- log2(data[GeneList, ]+1)
  MyCounts <- MyCounts - rowMeans(MyCounts)}
  else if (norm=="total")
  { NormData <- t(t(data)/total_counts(data)* mean(total_counts(data)))
  MyCounts <- log2(NormData[GeneList, ]+1)
  MyCounts <- MyCounts - rowMeans(MyCounts)}
  else if (norm=="DEseq")
  { NormData <- t(t(data)/median_of_ratios(data)) * mean(median_of_ratios(data))
  MyCounts <- log2(NormData[GeneList, ]+1)
  MyCounts <- MyCounts - rowMeans(MyCounts)}
  if (is.null(title))
  {
    title = "Heatmap based on log2 transformation"
  }
  pheatmap::pheatmap(MyCounts, main=title, cluster_rows=cluster_rows, cluster_cols=cluster_cols, annotation_col=annotation_col, annotation_row=annotation_row, fontsize_row= fontsize_row, fontsize_col= fontsize_col)
}

