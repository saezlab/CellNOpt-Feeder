\name{runDynamicFeeder}
\alias{runDynamicFeeder}

\title{
  Modelling the integrated PKN with CNORode
}

\description{
  This function evaluates the effects of possible feeder mechanisms which are added to the PKN.
}

\usage{
runDynamicFeeder(cnolist = cnolist, integratedModel = integratedModel, ode_parameters = ode_parameters, paramsSSm = paramsSSm, penFactor_k = 2)
}


\arguments{
  \item{cnolist}{a cnolist structure, as produced by \link{makeCNOlist}}
  \item{integratedModel}{the integrated model as returned from integrateLinks}
  \item{ode_parameters}{a list with the ODEs parameter information. Obtained with \link{createLBodeContPars}}
  \item{penFactor_k}{a penalty factor for the new integrated links obtained from the FEED algorithm and which are not present in the database (if the database was given). Default: penFactor_k = 100}
  \item{penFactor_tau}{a penalty factor for all the new nodes integrated in the PKN. Default: penFactor_tau = 1}
  \item{penFactorPIN_k}{a penalty factor for the new integrated links and which are present in the database (for the cases when the database was given). Default: penFactorPIN_k = 10}
  \item{paramsSSm}{a list of SSm parameters. default is the list returned by defaultParametersSSm}
}

\details{
This function evaluates the effects of possible feeder mechanisms which are added to the PKN. The analysis performed is a simple CNORode analysis over the integrated network where the newly integrated links are supposed to be penalised more than the links present in the original PKN. If a database of interactions is also provided by the user, than normally the links inferred from the FEED mechanism and which re not present in the database should be more penalised than the ones that are. There is also the opportunity to weight database interactions based on their relevance (i.e. number of resources, etc.). 
}

\value{
this function returns a list with fields:
\item{Parameters}{the inferred optimal ODE parameters}
\item{Integrated-Model}{the integrated model which was optimised}
\item{CNOlist}{the CNOlist object containing the data}

}

\author{
    E.Gjerga
}