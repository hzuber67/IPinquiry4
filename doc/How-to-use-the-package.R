## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(IPinquiry4)

## -----------------------------------------------------------------------------
CountTable <- system.file("extdata", "CountTable.txt", package = "IPinquiry4")
SampleTable <- system.file("extdata", "SampleTable.txt", package = "IPinquiry4")

## ----Function for data loading------------------------------------------------
IP_data <- load_IP_Data(CountTable, SampleTable)

## ----echo=TRUE, fig.height=5, fig.width=6-------------------------------------
MDSplot(IP_data)

## -----------------------------------------------------------------------------
#directory where is the file
IP_urt1 <- subset_IPObj_treat(IP_data, "urt1")

## ---- fig.height=5, fig.width=6-----------------------------------------------
MDSplot(IP_urt1, norm="DEseq")

## -----------------------------------------------------------------------------
test <- stat_test(IP_data, "urt1", treatment = "M1", div="DEseq", min.disp=10)

## -----------------------------------------------------------------------------
nrow(subset(test, test$adjp<0.05))

## -----------------------------------------------------------------------------
annotated_table_At <- addBiomaRtAnnotation(test, biomart="plants_mart", dataset="athaliana_eg_gene")

## ----results = 'hide'---------------------------------------------------------
createTable(annotated_table_At)

## -----------------------------------------------------------------------------
p <- createTable(annotated_table_At)
htmlwidgets::saveWidget(p,"interactive_table.html", selfcontained = TRUE)

## ---- fig.height=5, fig.width=7-----------------------------------------------
# add annotation as tag for the volcanoplot
htmlPlot(annotated_table_At, sign="adjp")

## ----results = 'hide'---------------------------------------------------------
# extract the 30 first letters of the description column
my_test = paste(row.names(annotated_table_At) , "_",substr(annotated_table_At$description,1,30))

# add annotation as tag for the volcanoplot
htmlPlot(annotated_table_At,  custom_text=my_test)

## ----results = 'hide'---------------------------------------------------------
smallTrueList <- c("AT1G26110.1", "AT5G45330.1", "AT3G13300")
htmlPlot(annotated_table_At, listGenes = smallTrueList, custom_text=my_test)

## ----results = 'hide'---------------------------------------------------------
# directory with the supplemental information for the example dataset
Supplemental <- system.file("extdata", "Supplemental_information.txt", package = "IPinquiry4")

# Add a supplemental column with criteria for color classification 
# based on the "Supplemental_information.txt" file
annotated_table_At2 <- add_suppl_information(annotated_table_At, Supplemental)

# Create the volcanoplot with colors according to this new column
htmlPlot(annotated_table_At2, colforcolor = annotated_table_At2$Classification)

## -----------------------------------------------------------------------------
p <- htmlPlot(annotated_table_At2, colforcolor = annotated_table_At2$Classification)
htmlwidgets::saveWidget(p,"interactive_volcanoplot_plot.html", selfcontained = TRUE)

## ----echo=TRUE, fig.height=5, fig.width=7-------------------------------------
PDF_Plot(annotated_table_At2)

## ----echo=TRUE, fig.height=5, fig.width=7-------------------------------------
PDF_Plot(annotated_table_At2, sign="p.value", max.pval = 0.05,  min.LFC = 1, line=TRUE, point_color= c("gray", "purple"), min_x=-10, max_x=10, min_y=0, max_y=20, point_size=3, label=TRUE, label_size=3, custom_text=annotated_table_At2$external_gene_name, title="The perfect volcanoplot")

## ----echo=TRUE, fig.height=5, fig.width=8-------------------------------------
PDF_Plot(annotated_table_At2, sign="adjp", max.pval = 0.05,  min.LFC = 1, line=TRUE, label=FALSE, colforcolor=annotated_table_At2$Classification, custom_text=annotated_table_At2$external_gene_name)

## ----echo=TRUE----------------------------------------------------------------
graph <- PDF_Plot(annotated_table_At2)
library(ggplot2)
ggsave("Volcanoplot.pdf", graph, height=7, width=10)

## ----echo=TRUE, fig.height=5, fig.width=8-------------------------------------
#Subset only protein with classification linked to RNA metabolism
class <- annotated_table_At2[!is.na(annotated_table_At2$Classification),]
# I create first a one-column dataframe with the classification for each selected proteins 
# and the protein ID as row.names
class2 <- class[,"Classification", drop=F]
# I draw the heatmap
IP_pheatmap(IP_data, GeneList=row.names(class), norm="DEseq", annotation_row = class2, fontsize_row=8)    

## -----------------------------------------------------------------------------
sessionInfo()

