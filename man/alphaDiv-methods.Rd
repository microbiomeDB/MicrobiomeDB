% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reexports-microbiomeComputations.R
\name{alphaDiv}
\alias{alphaDiv}
\title{Alpha diversity}
\usage{
alphaDiv(
  data,
  method = c("shannon", "simpson", "evenness"),
  verbose = c(TRUE, FALSE)
)
}
\arguments{
\item{data}{AbundanceData object}

\item{method}{string defining the the alpha diversity method. Accepted values are 'shannon','simpson', and 'evenness'}

\item{verbose}{boolean indicating if timed logging is desired}
}
\value{
ComputeResult object
}
\description{
This function returns alpha diversity values for each sample.
}
\examples{
alphaDivOutput <- alphaDiv(
     getCollection(
	microbiomeData::DiabImmune, 
	"16S (V4) Genus (Relative taxonomic abundance analysis)"), 
     method = "shannon"
)
}
