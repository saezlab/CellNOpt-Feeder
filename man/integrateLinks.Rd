\name{integrateLinks}
\alias{integrateLinks}

\title{
  Integrating the new links to the PKN
}

\description{
  This function integrates the new links inferred via the FEED method or from the database to the original PKN.
}

\usage{
integrateLinks(feederObject = feederObject, cnolist = cnolist, compression = FALSE, expansion = FALSE, database = database)
}


\arguments{
  \item{feederObject}{a feederObject structure, as produced by \link{buildFeederObjectDynamic}}
  \item{cnolist}{a cnolist structure, as produced by \link{makeCNOlist}}
  \item{compression}{a parameter indicating whether to compress the integrated model. Default: compression = FALSE}
  \item{expansion}{a parameter indicating whether to expand the integrated model. Default: expansion = FALSE}
  \item{database}{a database of interactions which can be optionally provided as an interaction matrix with 3 or 4 colums (source of interaction, sign of interaction, target of interaction and optionally a weight value from 0 to 1 indicating the significance of that interaction in the database). Default: database=NULL}
}

\details{
This function integrates the new links inferred via the FEED method or from the database to the original PKN. Moreover it indicates which are the integrated links and if a weighted database has been used it also shows the weights assigned to each integrated link. Links that are present in the original PKN are assigned a database weight of 0, integrated links that have been inferred via the FEED method and are not present in the database are assigned a database penalty of Inf, while integrated links present in the database take values between 0 and 1.
}

\value{
this function returns a list with fields:
\item{model}{the integrated model}
\item{integLinksIdx}{indices pointing towards the newly integrated links of the model}
\item{integSpeciesIdx}{indices pointing towards the newly integrated species of the model}
\item{databaseWeight}{weights assigned based on the presence of links in the database}

}

\author{
    E.Gjerga
}