\name{identifyMisfitIndices}
\alias{identifyMisfitIndices}

\title{
  Identifying indices of poorly fitted measurements
}

\description{
  This function identifies poorly fitted measurements for specific experimental conditions. It returns a list of possible indices and mse's pointing to possible connections to be added during the feeding process.
}

\usage{
identifyMisfitIndices(cnolist = cnolist, model = model, simData = simData, mseThresh = 0.05)
}


\arguments{
  \item{cnolist}{a cnolist structure, as produced by \link{makeCNOlist}}
  \item{model}{a model as returned by readSIF. Alternatively, the filename can also be provided.}
  \item{simData}{a matrix of simulated data values for a specific model as returned by plotLBodeFitness)}
  \item{mseThresh}{thrreshold parameter for minimal misfit to be considered (mseThresh = 0.05 by default)}
}
\details{
This function computes the misfits (MSE values) between the actual measured data points and the data values for a specific set of inferred model parameters. Once the MSE values are calculated for each of the measurements over each experimental condition, the poorly fitted measurements are then identify. A measurement is considered as poorly fitted if the corresponding inferred MSE value is higher than the specified MSE threshold value (mseThresh).
}

\value{
this function returns a list with fields:
\item{indices}{a list of indices pointing to the poorly fitted measurements and the corresponding ms value}
\item{use}{a matrix of use values indicating the mismatch between model simulations and data for each measurement at each experimental condition}
}
\author{
    E.Gjerga
}
\examples{
data(ToyModel_Gene, package="CNORfeeder")
data(CNOlistToy_Gene, package="CNORfeeder")
data(indices,package="CNORfeeder")
data(database, package="CNORfeeder")

indices = identifyMisfitIndices(cnolist = cnolist, model = model,
                                simData = simData, mseThresh = 0.05)
}