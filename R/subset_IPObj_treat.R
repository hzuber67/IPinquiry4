#'Subset IPinquiry object according to treatment
#'
#' Function to subset the IPObj according to the treatment column of the SampleTable
#'
#' @param IPObj object that comes directly from the function \code{load_IP_Data}
#' @param my_list List of the selected treatments (second column of the sampletable)
#' @return IPObj object for the selected treatments
#' @export

subset_IPObj_treat <- function(IPObj, my_list=c()) {
  ### Check arguments
  lapply(my_list, function(x) if (!x %in% levels(IPObj$treat)) { stop("Treatment name is not in the sampleTable") } )

  goodColumnsIDx <- c()
  for (i in my_list) {goodColumnsIDx <- c(goodColumnsIDx, which(IPObj$treat==i))}
  treat2 <- as.factor(as.character(IPObj$treat[goodColumnsIDx]))
  spectral_counts_2 <- IPObj$spectral_counts[,goodColumnsIDx]
  # Remove proteins not detected in any of the selected IPs
  spectral_counts_2 <- spectral_counts_2[apply(spectral_counts_2, 1, function(x) !all(x==0)),]
  if (length(IPObj) == 3) {batch2 <-as.factor(as.character(IPObj$batch[goodColumnsIDx]))
  C <- list(spectral_counts=spectral_counts_2 , treat=treat2, batch=batch2)
  C}
  else { C <- list(spectral_counts=spectral_counts_2 , treat=treat2)
  C}}
