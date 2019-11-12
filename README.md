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
install_github('saezlab/CellNOpt-Feeder', build_vignettes = TRUE)
# or download the source file from GitHub and install from source
install.packages('path_to_extracted_CellNOpt-Feeder-main_directory', repos = NULL, type="source")
```

## Toy Example
Here we show how Dynamic-Feeder can be applied on a single execution script over a simple toy example ([ToyMMB](http://www.bioconductor.org/packages/release/bioc/html/CellNOptR.html)).

Sourcing all the necessary packages
```R
library(CellNOptR)
library(MEIGOR)
library(CNORode)
library(doParallel)
library(readr)
library(infotheo)
library(igraph)
library(CNORfeeder)
```

Loading the toy example
```R
data("ToyModel", package="CellNOptR")
data("CNOlistToy", package="CellNOptR")

model = ToyModel
cnolist = CNOlist(CNOlistToy)
```

Setting the initial parameters (here parameters 'k' and 'tau' are optimised and 'n' fixed to 3) and optimization settings for the initial training of the model.
```R
ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0,
                                   LB_tau = 0, UB_n = 3, UB_k = 1, UB_tau = 1, default_n = 3,
                                   default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                   opt_tau = TRUE, random = TRUE)

## Parameter Optimization
# essm
paramsSSm=defaultParametersSSm()
paramsSSm$local_solver = "DHC"
paramsSSm$maxtime = 60;
paramsSSm$maxeval = Inf;
paramsSSm$atol=1e-6;
paramsSSm$reltol=1e-6;
paramsSSm$nan_fac=1000;
paramsSSm$dim_refset=30;
paramsSSm$n_diverse=1000;
paramsSSm$maxStepSize=Inf;
paramsSSm$maxNumSteps=10000;
paramsSSm$transfer_function = 4;

paramsSSm$lambda_tau=0.1
paramsSSm$lambda_k=0.01
paramsSSm$bootstrap=F
paramsSSm$SSpenalty_fac=0
paramsSSm$SScontrolPenalty_fac=0
```

Initial training of the model
```R
opt_pars=parEstimationLBode(cnolist, model, method="essm", ode_parameters=ode_parameters, paramsSSm=paramsSSm)
simData = plotLBodeFitness(cnolist = cnolist, model = model, ode_parameters = opt_pars, transfer_function = 4)
```

Identifying the mis-fits (measurements with mse worse than 0.05) and interactions from the database which we want to integrate (on this case only through data-driven method)
```R
indices = identifyMisfitIndices(cnolist = cnolist, model = model, simData = simData, mseThresh = 0.05)
feederObject = buildFeederObjectDynamic(model = model, cnolist = cnolist, indices = indices, database = NULL, DDN = TRUE)
integratedModel = integrateLinks(feederObject = feederObject, cnolist = cnolist)
```

Plotting the integrated model by highlighting in purple the new added links to the PKN
```R
plotModel(model = integratedModel$model, CNOlist = cnolist, indexIntegr = integratedModel$integLinksIdx)
```

Setting the initial ODE parameters to optimize for the integrated model
```R
ode_parameters=createLBodeContPars(integratedModel$model, LB_n = 1, LB_k = 0,
                                   LB_tau = 0, UB_n = 3, UB_k = 1, UB_tau = 1, default_n = 3,
                                   default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                   opt_tau = TRUE, random = TRUE)

```

Optimizing the integrated model with low penalty factor (lambda = 2) of the integrated links - here we observe the effects of the new links on the improvements in the fitting quality
```R
res1 = runDynamicFeeder(cnolist = cnolist, integratedModel = integratedModel, ode_parameters = ode_parameters, paramsSSm = paramsSSm, penFactor_k = 2)

plotLBodeFitness(cnolist = res1$CNOList, model = res1$`Integrated-Model`$model, ode_parameters = res1$Parameters, transfer_function = 4)
```

Optimizing the integrated model with high penalty factor (lambda_k = 10000) of the integrated links - here we do not observe any effects of the new links on the improvements in the fitting quality
```R
res2 = runDynamicFeeder(cnolist = cnolist, integratedModel = integratedModel, ode_parameters = ode_parameters, penFactor_k = 10000, paramsSSm = paramsSSm)

plotLBodeFitness(cnolist = res2$CNOList, model = res2$`Integrated-Model`$model, ode_parameters = res2$Parameters, transfer_function = 4)
```
