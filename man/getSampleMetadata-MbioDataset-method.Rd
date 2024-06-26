% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-MbioDataset.R
\name{getSampleMetadata,MbioDataset-method}
\alias{getSampleMetadata,MbioDataset-method}
\title{Get data.table of sample metadata from MbioDataset}
\usage{
\S4method{getSampleMetadata}{MbioDataset}(
  object,
  asCopy = c(TRUE, FALSE),
  includeIds = c(TRUE, FALSE),
  metadataVariables = NULL
)
}
\arguments{
\item{object}{MbioDataset}

\item{asCopy}{boolean indicating whether to return the data as a copy or by reference}

\item{includeIds}{boolean indicating whether we should include recordIdColumn and ancestorIdColumns}

\item{metadataVariables}{The metadata variables to include in the sample metadata. If NULL, all metadata variables will be included.}
}
\value{
data.table of sample metadata
}
\description{
Returns a data.table of sample metadata
}
\examples{
getSampleMetadata(microbiomeData::DiabImmune)
getSampleMetadata(microbiomeData::DiabImmune, metadataVariables = c("age_months", "sex"))
}
