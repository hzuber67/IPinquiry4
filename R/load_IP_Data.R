## Function for IP data loading
#' Load spectral count data and sample description.
#'
#' Function to load spectral count data. Spectral count for all IPs should be stored in one
#' tab delimited file with one distinct column for each IP experiments. The first column of the files should contains protein names.
#' The file needs to contain one header but must not contain any special characters. Name must not start by numbers.
#'
#' @param countFile tab delimited file (.txt) that contain all spectral counts with one distinct column for each IP experiments.
#' The first column of the files should contains protein names. The file needs to contain one header.
#' @param sampleTable tab delimited file (.txt) with two or three columns and an header. The first column correspond to the experiment names.
#' Each name should be unique and correspond to the names indicated in the header of \code{exprsFile}.
#' The second column indicates the different conditions (treatment, mutation). The third column is optional
#' and is used to indicate a batch effect.
#' @return \code{IPObj} object type that contains the dataframe with spectral counts and information about treatment and eventual batch effect
#' @export
#' @examples
#' CountTable <- system.file("extdata", "CountTable.txt", package = "IPinquiry4")
#' SampleTable <- system.file("extdata", "SampleTable.txt", package = "IPinquiry4")
#' IP_data <- load_IP_Data(CountTable, SampleTable)
#' @importFrom utils read.table

load_IP_Data <- function(countFile, sampleTable)
{
  #Create import table with spectral count
  Count_file <- read.table(countFile, sep= "\t", header = T, stringsAsFactors=default.stringsAsFactors(), quote="", fill=FALSE)
  Count_file_1 <- Count_file[,2:ncol(Count_file)]
  row.names(Count_file_1) <- Count_file[,1]
  Sample_table <- read.table(sampleTable, sep= "\t", header = T, stringsAsFactors=default.stringsAsFactors(), quote="", fill=FALSE)
  IPObj = list("spectral_counts"=Count_file_1)
  IPSampleTable <- Sample_table[,1]
  IPNames <- colnames(Count_file_1) # sample names of the count table
  IPCategSampleTable <- Sample_table[,2] # Conditions for each IP
  if (!is.unique(IPSampleTable))
  {stop("IP names are not unique. Change names of the \"first column\n",
        "of the SampleTable file")}
  # Check that number of samples are identical in the two files
  if (length(IPNames) != length(IPSampleTable))
  {stop("You must have the same number of samples in the SampleTable file \n",
        "as in the spectral count file")}
  # Tests to check that names on Count table and sample description table are indical. Order can be different between the two files
  if (!all(IPNames %in% IPSampleTable)) # On n'a pas les mêmes éléments dans les deux cas.
  {stop("Sample names must be the same in the SampleTable and the spectral count file")}
  # Check that there are at least two conditions
  if (length(unique(IPCategSampleTable)) == 1) {stop("You must have more than one condition to describe your experiments")}
  # Warning in case there are more than two conditions
  else if (length(unique(IPCategSampleTable)) > 2){warning("You have more than two conditions. Statistical analysis is possible with only two factors. You have to fill the \"treatment\" option")}
  # Make sample order similar
  IPObj$treat <- IPCategSampleTable[sapply(IPNames, function(x) which(IPSampleTable==x))]
  #Il peut y avoir une troisième colonne, qui sera celle du batch effect.
  if (length(colnames(Sample_table))==3)
  {
    IPBatchSampleTable <- Sample_table[,3]
    IPObj$batch <- IPBatchSampleTable[sapply(IPNames, function(x) which(IPSampleTable==x))]
  }
  IPObj
}
