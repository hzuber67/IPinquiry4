# IPinquiry4
## Statistical analysis and data visualisation for co-immunoprecipitation coupled to mass spectrometry experiments

**IPinquiry** is a package written for **R** to analyze and visualize data from co-immunoprecipitation (co-IP) coupled to mass spectrometry experiments. In particular, IPinquiry compiles several functions to:

- Identify protein significanlty enriched in co-IP experiments
- Retrieve annotation for detected proteins
- Export result tables
- Create different graph types, such as volcano plot, multidimensional scaling or heatmap.

The main purpose of this package is to provide a simple R pipeline with a limited number of processing steps to facilitate as much possible data exploration and plot creation for biologists.

For the statistical analysis, IPinquiry uses the **Genewise Negative Binomial Generalized Linear Model** developped by [**EdgeR**](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and uses **spectral counts** as input. This statistical model is also among the models proposed in [**msmsTests**](https://www.bioconductor.org/packages/release/bioc/html/msmsTests.html) to analyze mass spectrometry data .  


### Package requirements
To use **IPinquiry4**, you first need to download and install R. We also recommend the use of [**RStudio**](https://rstudio.com), which provides a nice R user interface making life easier for R beginners. Also, [**EdgeR**](https://bioconductor.org/packages/release/bioc/html/edgeR.html) is required for IPinquiry4 installation and the following R packages are needed for one or several IPinquiry4 functions :

- ggplot2,
- xlsx,
- RColorBrewer,
- pheatmap,
- DT,
- plotly,
- biomaRt,
- htmlwidgets


### Install IPinquiry4

IPinquiry4 can be downloaded and installed from Github using [**devtools**](https://cran.r-project.org/web/packages/devtools/index.html).

Install devtools (if needed) and load the library
```{R}
install.packages("devtools")
library(devtools)
```

IPinquiry4 can then be installed. Add the argument `build_vignettes = TRUE`, if you want to be able acess the html vignette "How to use the IPinquiry4 R package".

```{R}
install_github("https://github.com/hzuber67/IPinquiry4", build_vignettes = TRUE)
```


### Run IPinquiry4

Please refer to the html tutorial for a detailed IPinquiry protocol. An example dataset is included in the packages.

```{R}
browseVignettes("IPinquiry4")
```
Also, documentation for each individual function can be loaded using `help('the function')`
