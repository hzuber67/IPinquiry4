#' Subset IPinquiry object according to batch
#'
#' Function to subset the IPObj according to the batch column of the SampleTable
#'
#' @param IPObj IPObj type object that comes directly from the function \code{load_IP_Data}
#' @param my_list List of the selected batchs (third column of the sampletable)
#' @return IPObj object with only the selected batchs
#' @export
subset_IPObj_batch <- function(IPObj, my_list=c()) {
  ### Check arguments
  lapply(my_list, function(x) if (!x %in% levels(IPObj$batch)) { stop("Batch name is not in the sampleTable") } )
  goodColumnsIDx <- c()
  for (i in my_list) {goodColumnsIDx <- c(goodColumnsIDx, which(IPObj$batch==i))}
  batch2 <- as.factor(as.character(IPObj$batch[goodColumnsIDx]))
  spectral_counts_2 <- IPObj$spectral_counts[,goodColumnsIDx]
  treat2 <- as.factor(as.character(IPObj$treat[goodColumnsIDx]))
  spectral_counts_2 <- spectral_counts_2[apply(spectral_counts_2, 1, function(x) !all(x==0)),]
  C <- list(spectral_counts=spectral_counts_2 , treat=treat2, batch=batch2)
  C}

