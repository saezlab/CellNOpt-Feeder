# CellNOpt-Feeder

CNORfeeder is an add-on to [CellNOptR](www.cellnopt.org) that permits to extend
a network derived from literature with links derived in a strictly data-driven 
way and supported by protein-protein interactions as described in:

F. Eduati, J. De Las Rivas, B. Di Camillo, G. Toffolo, J. Saez-Rodriguez. [Integrating literature-constrained and data-driven inference of signalling networks](http://bioinformatics.oxfordjournals.org/content/28/18/2311). *Bioinformatics*, 2012 , 28 (18)

## Installation

CNORfeeder is currently available for the installation as an R-package from our GitHub page

```R
# Install CARNIVAL from Github using devtools
# install.packages('devtools') # in case devtools hasn't been installed
library(devtools)
install_github('saezlab/CellNOpt-Feeder')
# or download the source file from GitHub and install from source
install.packages('path_to_extracted_CellNOpt-Feeder_directory', repos = NULL, type="source")
```
