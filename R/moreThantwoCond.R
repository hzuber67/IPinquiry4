#' Create sub dataset
#'
#' \code{moreThantwoCond} function called by the statistical test function (\code{modifFINALEdgeR})
#' to reduce the dataset to two conditions when more are in the input table.
#' @param IPObj Object loaded by the loadData function
#' @param cond1 Name of the first condition (control)
#' @param cond2 Name of the second condition
#' @param batchE Determine if a batch effect in taken in the statistical model
#' @return Return an IPObj type object that contains only the count for the selected conditions
#' @importFrom utils write.table
#'
moreThantwoCond <- function(IPObj, cond1, cond2, batchE=FALSE)
{
  ## create the argument facs= sample table
  if (batchE)
  {facs <- as.matrix(data.frame(IPObj$treat, IPObj$batch))
  colnames(facs) <- c("treat","batch")
  row.names(facs) <- colnames(IPObj$spectral_counts)}
  else
  {facs <- as.matrix(data.frame(IPObj$treat))
  colnames(facs) <- "treat"
  row.names(facs) <- colnames(IPObj$spectral_counts)}

  dat <- as.data.frame(IPObj$spectral_counts)
  goodColumnsIDx1 <- sapply(cond1, function(x) which(IPObj$treat==x))
  goodColumnsIDx2 <- sapply(cond2, function(x) which(IPObj$treat==x))
  goodColumnsIDx <- c(goodColumnsIDx1, goodColumnsIDx2)
  # Create a temporary subtable
  ncols = 2
  tmpVec <- c(rownames(facs)[goodColumnsIDx], as.vector(IPObj$treat[goodColumnsIDx]))
  if (length(colnames(facs))==2)# il y a batch effect
  {
    ncols=3
    tmpVec <- c(tmpVec,as.vector(IPObj$batch[goodColumnsIDx]))
  }
  tmpsample <- as.data.frame(matrix(tmpVec,ncol=ncols),row.names=NULL)
  allIDx <- 1:length(colnames(dat))
  badIDs <- allIDx[sapply(allIDx, function(x) !x %in% goodColumnsIDx)]
  dat[,badIDs] <- NULL
  dat$names <- rownames(dat)
  lenDat <- length(colnames(dat))
  ahhh <-dat[,c(lenDat, 2:lenDat-1)]
  nameTmpFile <- paste(MHmakeRandomString(), ".txt", sep="")
  nameTmpsample <- paste(MHmakeRandomString(), ".txt", sep="")
  write.table(tmpsample, nameTmpsample , sep = "\t", col.names = T, row.names=F, quote=F)
  write.table(dat[,c(lenDat, 2:lenDat-1)], nameTmpFile, sep = "\t", col.names = T, row.names=F, quote=F)
  newMSN <- load_IP_Data(nameTmpFile,nameTmpsample)
  file.remove(nameTmpFile)
  file.remove(nameTmpsample)
  newMSN
}
