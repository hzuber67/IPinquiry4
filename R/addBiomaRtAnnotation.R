#' Automatic annotation using BiomartR package
#'
#' \code{addBiomaRtAnnotation} is a function to retrieve protein, transcript or gene annotation
#'
#' @param normObj dataframe with row.names that contains gene, protein or transcript ID
#' @param biomart \code{plants_mart} by default (for plants). For other organisms, use  "ENSEMBL_MART_ENSEMBL",
#' see \code{listMarts} biomartR function
#' @param dataset \code{"athaliana_eg_gene"} by default. biomartR dataset that corresponds to the studied organism.
#' List of available datasets can be obtained using  biomartR \code{listDatasets} function
#' @param host \code{"www.plants.ensembl.org"} by default. Can be "www.ensembl.org" to get other dataset.
#' @param features \code{"ensembl_peptide_id"}, \code{"ensembl_transcript_id"}, \code{"ensembl_gene_id"}, \code{"external_gene_name"}
#' \code{"ensembl_peptide_id"} by default. To be adjusted according to ID in the row.names of the table
#' @return Return a result table with additional columns that contain gene or protein annotation.
#' @export
addBiomaRtAnnotation_2 <- function(normObj, biomart="plants_mart",
                                   dataset="athaliana_eg_gene", host="plants.ensembl.org",
                                   features=c("ensembl_peptide_id","ensembl_transcript_id","ensembl_gene_id", "external_gene_name"))
{
  stopifnot(!missing(normObj))
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package \"biomaRt\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  filter <- match.arg(features)
  mart <- biomaRt::useMart(biomart,dataset,host=host)
  if (features=="external_gene_name") {
    atr <- c(filter,"ensembl_gene_id","description")
  }
  else {
    atr <- c(filter,"external_gene_name","description")
  }
  geneInfo <- biomaRt::getBM(rownames(normObj),mart = mart, attributes = atr,
                             filters = filter) # take only ID that match feature found on ensembl.org
  if (length(geneInfo$description)==0) {stop(sprintf("rowIDs from count table are different than IDs from \"%s\" from ensembl\nChoose another feature ID, or another dataset from biomaRt, or check that your IDs are in the format as the one from ensembl", filter))}
  if (length(geneInfo$description) < length(normObj[,1])){
    warning("Annotation was not found for all ID\n")}
  if (length(duplicated(geneInfo))>0){
    warning("Duplicate annotations were found for some ID => these annotations are merged\n")
    geneInfo <- aggregate(.~geneInfo[,filter], data=geneInfo, paste, collapse = "|")
  }
  rownames(geneInfo) <- geneInfo[,1] # first column of geneInfo table contain gene ID.
  geneInfo[,1] <- NULL
  normObj <- merge(normObj, geneInfo, by="row.names", all.x=TRUE)
  rownames(normObj) <- normObj$Row.names
  normObj$Row.names <- NULL
  if ("number" %in% colnames(normObj)) {normObj <- normObj[order(normObj$number),] }
  normObj
}
