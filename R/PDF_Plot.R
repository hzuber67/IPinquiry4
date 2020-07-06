################################################################################
##               Create volcanoplot graph suitable for pdf export                      #
################################################################################
#' Create volcanoplot graph suitable for pdf export
#'
#' \code{PDF_Plot} create a volcano plot using the ggplot2 package
#'
#' @param tres table obtained from the  \code{stat_test} function, that is one table with one ID by row and
#' at least one column with \code{logFC} and one other with  \code{adjp} or \code{p.value}
#' @param sign Choice a the value used to show statistical significance: \code{adjp} (default) is \code{p.value}
#' @param max.pval Cut-off p-value considered as relevant
#' @param min.LFC Minimum LogFC for a protein to be considered as significant
#' @param line Dotted red lines to show FC and p-value cut-off. it can be \code{TRUE} or \code{FALSE}
#' @param min_x Optional, set min value for x axis
#' @param max_x Optional, set max value for x axis
#' @param min_y Optional, set min value for y axis
#' @param max_y Optional, set max value for y axis
#' @param colforcolor Optional, colored protein according to information in a column f the table
#' @param custom_text Optional, allow for setting custom text in labels
#' @param label Set \code{FALSE} to remove label
#' @param label_size Optional, allow for manually setting label sizes
#' @param point_color Optional, allow for manually setting colors
#' @param point_size Optional, allow for manually setting point sizes
#' @param title Optional, allow for adding a title
#' @return Return graph that is suitable for pdf export
#' @details
#' If \code{colforcolor} is empty, protein with significantly different FC will be highlighted
#' If \code{colforcolor}, is filled, point will be colored according to information in the indicated column
#' @export

PDF_Plot <- function (tres, sign=c("adjp", "p.value"), max.pval = 0.05,  min.LFC = 1, line=TRUE, point_color= NULL, min_x=NULL, max_x=NULL, min_y=NULL, max_y=NULL, colforcolor=NULL, point_size=4, label=TRUE, label_size=4, custom_text=NULL, title=NULL)
{

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
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


  if (is.null(point_color)) {

    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("Package \"RColorBrewer\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    point_color=RColorBrewer::brewer.pal(12, "Paired")

    }

  if (is.null(custom_text)) {tres$names <- row.names(tres)}
  else if (!is.null(custom_text)) {tres$names <- custom_text}

  #definition y et x maximum
  if (is.null(min_y) || is.null(max_y)) {
    min_y = 0
    max_y = max(-log10(tres$sign)) + 0.1 * (max(-log10(tres$sign)))
  }

  if (is.null(min_x) || is.null(max_x)) {
     min_x = min(tres$LogFC) - 0.1 * (max(tres$LogFC) - min(tres$LogFC))
     max_x = max(tres$LogFC) + 0.1 * (max(tres$LogFC) - min(tres$LogFC))
  }
  #make plot
  LogFC <- filter <- NULL
  volcano <- ggplot2::ggplot(tres, ggplot2::aes(x=LogFC, y=-log10(sign)))

  #graph par defaut color according to p.value
  if(is.null(colforcolor))
  {
    if ("filter" %in% colnames(tres))
    {
      volcano <- volcano +
        ggplot2::geom_point(ggplot2::aes(colour=sign<max.pval&filter == "YES"), size=point_size, shape=19) +
        ggplot2::scale_color_manual(name = 'significant p-value', values=point_color)
    }
    else
    {
      volcano <- volcano +
        ggplot2::geom_point(ggplot2::aes(colour=sign<max.pval), size=point_size, shape=19) +
        ggplot2::scale_color_manual(name = 'significant p-value', values=point_color)
    }
  }

  #graph with color according to other criteria
  else if(!is.null(colforcolor))
  {
    if (sum(is.na(colforcolor)) > 0 && is.character(colforcolor))
    { #colforcolor <- as.character(colforcolor)
      colforcolor[is.na(colforcolor)]<- "no_info" }
    else if (sum(is.na(colforcolor)) > 0 && is.numeric(colforcolor))
    { #colforcolor <- as.nueric(colforcolor)
      colforcolor[is.na(colforcolor)]<- 0 }

    volcano <- volcano + ggplot2::geom_point(ggplot2::aes(color=colforcolor), size=point_size, shape=19) + ggplot2::scale_color_manual(name="information", values=point_color)

  }

  #add label for significant protein

  if ("filter" %in% colnames(tres) && label)
  {
    sub <- tres[tres$sign<0.05&tres$filter == "YES",]
    volcano <- volcano +
      ggplot2::geom_text(data=sub, label=sub$names, vjust = 0, nudge_y = 0.2, size= label_size)

  }
  else if (label)
  {
    sub <- tres[tres$sign<0.05,]
    volcano <- volcano +
      ggplot2::geom_text(data=sub, label=sub$names, vjust = 0, nudge_y = 0.2, size= label_size)

  }
  # Add dotted red lines to indicate FC and p-value cut-offs
  if (line) {

    volcano <- volcano + ggplot2::geom_hline(yintercept = -log10(max.pval), linetype="dashed", color = "red") +
                         ggplot2::geom_vline(xintercept = min.LFC, linetype="dashed", color = "red") +
                         ggplot2::geom_vline(xintercept = -min.LFC, linetype="dashed", color = "red")
  }

  #Make final format of the plot
  volcano <- volcano + ggplot2:: scale_x_continuous(limits=c(min_x , max_x)) + ggplot2::scale_y_continuous(limits=c(min_y , max_y)) +
                       ggplot2::theme(
                                      axis.line.x = ggplot2::element_line(colour = "black", size=0.3),
                                      axis.line.y = ggplot2::element_blank(),
                                      axis.text=ggplot2::element_text(size=12),
                                      axis.ticks=ggplot2::element_line(size=0.8),
                                      panel.grid.major = ggplot2::element_blank(),
                                      panel.grid.minor = ggplot2::element_blank(),
                                      panel.background = ggplot2::element_blank(),
                                      legend.text = ggplot2::element_text(size=10),
                                      legend.title = ggplot2::element_text(size=12),
                                      axis.title = ggplot2::element_text(size=12),
                                      plot.title = ggplot2::element_text(hjust = 0.5, face="bold")
                                      ) +
                       ggplot2::xlab("LogFC") +
                       ggplot2::ylab(cust_title) +
                       ggplot2::ggtitle (title) +
                       ggplot2::geom_vline(xintercept = 0, color = "black", size=0.3)
  volcano
}

