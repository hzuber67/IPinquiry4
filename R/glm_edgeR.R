#' Statistical model from the EdgeR package

#' \code{batchVerification} Function to check that batch was indicated in the SampleTable
#' @param IPObj object, that comes from the function \code{loadData}
#' @return Return \code{TRUE} if column with batch information exists
batchVerification <- function(IPObj)
{
  if (is.null(IPObj$batch)) {stop("IPObj object does not have a batch column")}
  TRUE
}

#' Genewise negative binomial generalized linear model developped by the EdgeR
#'
#' \code{modifFINALEdgeR} Function
#'
#' @param IPObj object from function \code{loadData}.
#' @param control Name of the control. The name must be concordant with one of the conditions' name of the sample
#' table (2nd column)
#' @param treatment \code{NULL} by default, will take the name of the second condition in the level order. Can be
#' changed by a name included in the conditions column of the sample table
#' @param batchE \code{FALSE} by default. If there is batch effect in the third column of the sample table, can be set to \code{TRUE}
#' @param fnm \code{NULL} or a character string with the treatment factor name,
#' as used in the column names of the factors data frame, and in the formula.
#' @param div Scale factor
#' @param filter Minimal number of spectra
#' @param abundance.trend argument of EdgeR function
#' @param min.disp Minimal number to calculate dispersion
#' @param glm Genewise Negative Binomial Generalized Linear Model with \code{"QL"} or without \code{"classic"} Quasi-likelihood tests.
#' @return Return tables with four column
#' @importFrom stats relevel model.matrix as.formula

glm_edgeR <- function (IPObj,control,treatment=NULL,batchE=FALSE, div=NULL,fnm=NULL, filter=NULL, min.disp=NULL, abundance.trend=TRUE, glm=NULL)
{
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package \"edgeR\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if(!is.null(treatment))
  {
    if(!treatment %in% IPObj$treat){stop("Error in treatment name :\ntreatment must be one of the condition written in sampleTable")}
    if(!is.null(treatment) && length(unique(IPObj$treat)) > 2 && (batchE) ){IPObj <- moreThantwoCond(IPObj, control, treatment, batchE)}
    if (!is.null(treatment) && length(unique(IPObj$treat)) > 2 && (!batchE)){IPObj <- moreThantwoCond(IPObj, control, treatment)}
  }
  else if (length(unique(IPObj$treat)) > 2) {stop("You have more than two conditions, you have to fill the \"traitement\" argument")}

  data <- data.matrix(IPObj$spectral_counts)
  data <- data[apply(data, 1, function(x) !all(x==0)),]

  #normalisation factor calculation
  if (div== "nothing") {div <- NULL}
  else if (div== "total") {div <- total_counts(IPObj$spectral_counts)}
  else if (div== "DEseq") {div <- median_of_ratios(IPObj$spectral_counts)}

  #numeric scalar giving a value for the filtering out of low abundance tags.
  if (is.null(min.disp)) {min.disp <- 5} # defaut value proposed by EdgeR
  ## create the argument facs= sample table
  if (batchE)
  {facs <- as.matrix(data.frame(IPObj$treat, IPObj$batch))
  colnames(facs) <- c("treat","batch")
  row.names(facs) <- colnames(IPObj$spectral_counts)}
  else
  {facs <- as.matrix(data.frame(IPObj$treat))
  colnames(facs) <- c("treat")
  row.names(facs) <- colnames(IPObj$spectral_counts)}

  if (is.null(div))
    div <- rep(1, ncol(data))
  if (is.null(fnm))
    fnm <- colnames(facs)[1]
  if (is.null(IPObj$treat)) {stop("\"treat\" variable from IPObj object is empty")}
  if (!control %in% IPObj$treat) {stop("Declared control reference is not in the sampleTable")}
  if (is.null(treatment) && length(unique(levels(IPObj$treat)))>2) {stop("You have to fill the \"treatment\" parameter because you have more than two conditions")}
  tmpDataFrame <- data.frame(facs)
  if (!is.null(treatment))
  {
    if (!treatment %in% IPObj$treat) {stop("The second declared condition is not in sampleTable")}
    else {tmpDataFrame$treat <- relevel(as.factor(tmpDataFrame$treat), treatment)}
  }
  tmpDataFrame$treat <- relevel(as.factor(tmpDataFrame$treat), control)
  if (batchE) # We have a batch effect
  {
    if(batchVerification(IPObj))#,facs))
    {
      M <- model.matrix(~tmpDataFrame$treat+tmpDataFrame$batch)
      Mo <- model.matrix(~tmpDataFrame$batch)
    }
  }
  else # Don't have a batch effect
  {
    M <- model.matrix(~tmpDataFrame$treat)
    Mo <- model.matrix(as.formula("~1"), data.frame(facs))
  }
  y <- edgeR::DGEList(counts = data, group = facs[, fnm])
  y$offset <- log(div)
  y <- edgeR::estimateDisp(y, min.row.sum=min.disp, design = M)
  y$offset <- log(div)

  if (is.null(glm)) {glm <- "QL"}
  if( glm=="QL")
  {
    fit <- edgeR::glmQLFit(y, design = M, abundance.trend=abundance.trend, robust=TRUE)
    vc <- setdiff(colnames(M), colnames(Mo))
    res <- edgeR::glmQLFTest(fit, coef = vc)
  }
  else if (glm=="classic")
  {
    fit <- edgeR::glmFit(y, design = M)
    vc <- setdiff(colnames(M), colnames(Mo))
    res <- edgeR::glmLRT(fit, coef = vc)
  }

  tres <- res$table[, c(1, 3, 4)]
  if( glm=="QL")
  {colnames(tres) <- c("LogFC", "F.statistic", "p.value")}
  if( glm=="classic")
  {colnames(tres) <- c("LogFC", "LR", "p.value")}
  if (!is.null(filter))
  {
    print("Supplemental filter")
    tres$filter <- ifelse(rowSums(data > filter) >= 1,"YES", "NO")
  }
  tres
}

