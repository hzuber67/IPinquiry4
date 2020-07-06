#' Statistical analysis based on the glm model developped by EdgeR
#'
#' Function to perform the statistical analysis using the genewise negative binomial generalized linear model developped by the EdgeR package
#' By defaut, IPinquiry4 package uses GLM model with Quasi-likelihood Tests. By adding the argument `glm="classic"`,
#' you can use instead the EdgeR function based on the GLM without Quasi-likelihood Tests.

#' The \code{stat_test} function compares data from two conditions.
#' @param IPObj object from the function \code{load_IP_Data}.
#' @param control denominator, name of the condition used as a reference. This name should correspond to one of the condition of the second
#' column of the \code{sampleTable}.
#' @param treatment Numerator, name of the condition to be compared with the \code{control}.
#' By default, the second condition found in the \code{sampleTable} is used as treatment.
#' @param batch Boolean-type parameter (default FALSE), if TRUE and if the \code{sampleTable} contains a third table
#' with batch indications, the batch effect are taken into account in the statistical model.
#' @param div Choice of the method of calculation for the scale factor
#' @param min.disp numeric scalar giving a value for the filtering out of low abundance tags for the calculation of the dispersion.
#' Only tags with total sum of counts above this value are used. Low abundance tags can adversely affect the dispersion estimation,
#' so this argument allows the user to select an appropriate filter threshold for the tag abundance. By defaut, this value is set to 5 by EdgeR.
#' @param filter Supplemental filter. Minimal number of spectra to be considered as significant, indicated as a supplemental column with YES or NO
#' @param glm Genewise Negative Binomial Generalized Linear Model with \code{"QL"} or without \code{"classic"} Quasi-likelihood tests.
#' @return Return a dataframe with five colmumns: \code{logFC} (control used as denominator), \code{LR},
#' \code{p.value}, \code{adjp} and \code{number} corresponding the rank according to ajusted p-value.
#' The adjusted p-value is calculated according to the False Discovery Rate (FDR) method.
#' @details If \code{div} is \code{NULL} or \code{nothing}, scales factors are set to 1.
#' If \code{div} is \code{total}, scales factors are the total number of counts by sample.
#' If \code{div} is \code{DEseq}, scales factors are calculated according to the "median-to-ratio" method as in DEseq.
#' @examples
#' CountTable <- system.file("extdata", "CountTable.txt", package = "IPinquiry4")
#' SampleTable <- system.file("extdata", "SampleTable.txt", package = "IPinquiry4")
#' IP_data <- load_IP_Data(CountTable, SampleTable)
#' test <- stat_test(IP_data, "urt1", treatment = "M1", div="DEseq")
#' @export
#' @importFrom stats p.adjust

stat_test <- function(IPObj, control, treatment=NULL, batch=FALSE, div=c("nothing", "total", "DEseq"), filter=NULL, min.disp=NULL, glm=NULL)
{
  ### I Check argument
  div <- match.arg(div)

  if (!control %in% levels(IPObj$treat))
  {
    stop("Error in control name :\ncontrol must be one of the condition written in sampleTable")
  }

  if (!is.null(treatment))
  {
    if (!treatment %in% levels(IPObj$treat))
    {stop("Error in treatment name :\ntreatment must be one of the condition written in sampleTable")}
    ## d : control et traitement doivent être différents.
    if (control == treatment)
    {stop("control and treatment must be different")}
  }

  if (batch)
  {
    if(!"batch" %in% names(IPObj))
    {
      stop("You have to declare the batch effect in the sampleTable, third column")
    }
  }
  if (is.null(glm)) {glm <- "QL"}
  ### II Call the EdgeR glm function and adjust p-value
  normObj <- glm_edgeR(IPObj,control,treatment=treatment, batchE=batch, div=div, filter=filter, min.disp=min.disp, abundance.trend=TRUE, glm=glm)
  normObj$adjp <- p.adjust(normObj$p.value, method = "fdr")
  normObj <- normObj[order(normObj$adjp),]
  normObj$number <- c(1:nrow(normObj))
  normObj
}
