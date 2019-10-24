# CNORfeeder

CNORfeeder is an add-on to [CellNOptR](www.cellnopt.org) that permits to extend
a network derived from literature with links derived in a strictly data-driven 
way and supported by protein-protein interactions as described in:

F. Eduati, J. De Las Rivas, B. Di Camillo, G. Toffolo, J. Saez-Rodriguez. [Integrating literature-constrained and data-driven inference of signalling networks](http://bioinformatics.oxfordjournals.org/content/28/18/2311). *Bioinformatics*, 2012 , 28 (18)

## Requirements

All network modelling steps were performed in [R](https://www.r-project.org/) v3.5.1 and visualised using [Cytoscape v3.4.](http://chianti.ucsd.edu/cytoscape-3.4.0/)

Other prerequisites include downloading and installing the following `R` package dependencies:

### Bioconductor

+ [CellNOptR](https://bioconductor.org/packages/release/bioc/html/CellNOptR.html)
+ [MEIGOR](https://www.bioconductor.org/packages/release/bioc/html/MEIGOR.html)
+ [CNORode](https://github.com/saezlab/CNORode)

### CRAN
+ [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)
+ [readr](https://cran.r-project.org/web/packages/readr/index.html)
+ [infotheo](https://cran.r-project.org/web/packages/infotheo/index.html)
+ [igraph](https://cran.r-project.org/web/packages/igraph/index.html)

## Installation
CellNOpt-Feeder is currently available for the installation as an R-package from our GitHub page

```R
# Install CARNIVAL from Github using devtools
# install.packages('devtools') # in case devtools hasn't been installed
library(devtools)
install_github('saezlab/CellNOpt-Feeder')
# or download the source file from GitHub and install from source
install.packages('path_to_extracted_CellNOpt-Feeder-main_directory', repos = NULL, type="source")
```
