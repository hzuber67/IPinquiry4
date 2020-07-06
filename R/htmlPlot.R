#' Interactive volcanoplot using the R package Plotly
#'
#' \code{htmlPlot} creates an interactive volcanoplot showing statistical significance (p-value or adjusted p-value) vs protein enrichment (LogFC)
#'
#' @param tres Table obtained from the \code{stat_test} function, that is one table with one ID by row and
#' at least one column with \code{logFC} and one other with  \code{adjp}
#' @param sign Choice a the value used to show statistical significance: \code{adjp} (default) is \code{p.value}
#' @param max.pval Cut-off of p-value considered as relevant
#' @param min.LFC Minimum LogFC for a protein to be considered as significant
#' @param listGenes Optional, list of gene/protein ID to be colored
#' @param colforcolor Optional, colored protein according to information in a supplemental table
#' @param custom_text Optional, add text labels
#' @return Return html interactive graph
#' @details
#' If \code{listGenes} is empty, protein with significantly different FC will be highlighted
#' If \code{listGenes}, is filled, only proteins contained in the list will be filled
#' @export

htmlPlot <- function (tres, sign=c("adjp", "p.value"), max.pval = 0.05, min.LFC = 1, listGenes = c(), colforcolor=NULL, custom_text=NULL)
{
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package \"plotly\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  sign <- match.arg(sign)
  if (sign=="adjp") {

    tres$sign <- tres$adjp
    cust_title <- "-log10(adjp)"
  }

  if (sign=="p.value") {

    tres$sign <- tres$p.value
    cust_title <- "-log10(p.value)"
  }

  if (is.null(custom_text)) {tres$names <- row.names(tres)}
  else if (!is.null(custom_text)) {tres$names <- custom_text}
  min_y = 0
  max_y = max(-log10(tres$sign)) + 0.1 * (max(-log10(tres$sign)))
  min_x = min(tres$LogFC) - 0.1 * (max(tres$LogFC) - min(tres$LogFC))
  max_x = max(tres$LogFC) + 0.1 * (max(tres$LogFC) - min(tres$LogFC))
  if ("filter" %in% colnames(tres))
  { pflags <- tres$sign <= max.pval & abs(tres$LogFC) >= min.LFC & tres$filter == "YES"}
  else
  { pflags <- tres$sign <= max.pval & abs(tres$LogFC) >= min.LFC}
  p <- plotly::plot_ly()# %>%
  if(length(listGenes)>0 && !is.null(colforcolor))
  {
    stop("You have to chose between /listGenes and /colforcolor arguments")
  }
  if(length(listGenes)==0 && is.null(colforcolor))
  {
    p<-plotly::add_trace(p,tres, x = tres$LogFC[!pflags], y = -log10(tres$sign[!pflags]), mode = "markers", type="scatter", marker = list(size=10), text = tres$names[!pflags], name = "no_signif")
    p<-plotly::add_trace(p,tres, x = tres$LogFC[pflags], y = -log10(tres$sign[pflags]), mode = "markers", type="scatter", marker = list(size=10), text = tres$names[pflags], name = "signif")
  }
  else if(!is.null(colforcolor))
  {
    if (sum(is.na(colforcolor)) > 0 && is.character(colforcolor))
    { #colforcolor <- as.character(colforcolor)
      colforcolor[is.na(colforcolor)]<- "no_info" }
    else if (sum(is.na(colforcolor)) > 0 && is.numeric(colforcolor))
    { #colforcolor <- as.nueric(colforcolor)
      colforcolor[is.na(colforcolor)]<- 0 }

    p<-plotly::add_trace(p,tres, x = tres$LogFC, y = -log10(tres$sign), color = colforcolor, type="scatter", mode = "markers", marker = list(size=10), text = tres$names, name = colforcolor)
    p<-plotly::add_trace(p,tres, x = tres$LogFC[pflags], y = -log10(tres$sign[pflags]), mode = "markers", type="scatter", marker = list(size=10, color='rgba(17, 157, 255,0)',line = list(color = 'black',width = 1.5)), text = tres$names[pflags], name = "signif")
  }
  else # Une liste avec des identifiants d'interêt
  {
    interestData <- tres[c(listGenes),]
    if ("NA" %in% rownames(interestData))
    {warning("One or severall ID of your interst list are not in the original dataframe")}
    p<-plotly::add_trace(p,tres, x = tres$LogFC, y = -log10(tres$sign), type="scatter", mode = "markers", marker = list(size=10), text = tres$names, name = "no_interest")
    p<-plotly::add_trace(p,interestData, x = interestData$LogFC, y = -log10(interestData$sign), type="scatter", mode = "markers", marker = list(size=10), text = interestData$names, name = "interset")
  }
  p<-plotly::add_trace(p,x=c(min_x, max_x), y=c(-log10(max.pval),-log10(max.pval)) , type="scatter", mode="lines", name='hline', line = list(dash= "dash", color = "red"))#%>% # ligne horizonale
  p<-plotly::add_trace(p,x=c(-min.LFC, -min.LFC), y = c(0, max_y), type="scatter", mode="lines", name='vline1', line = list(dash= "dash", color = "red"))#%>%
  p<-plotly::add_trace(p,x=c(min.LFC, min.LFC), y = c(0, max_y), type="scatter", mode="lines", name='vline2', line = list(dash= "dash", color = "red"))#%>%
  p<-plotly::layout(p,title = 'Volcano plot',
            yaxis = list(
              title = cust_title,
              range = c(min_y,max_y)),
            xaxis = list(
              title = "logFC",
              range = c(min_x, max_x))
  )
  p
}
