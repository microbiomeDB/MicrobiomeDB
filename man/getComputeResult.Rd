% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{getComputeResult}
\alias{getComputeResult}
\alias{getComputeResult,ComputeResult-method}
\alias{getComputeResult,CorrelationResult-method}
\alias{getComputeResult,DifferentialAbundanceResult-method}
\title{Get Microbiome Dataset Compute Result}
\usage{
getComputeResult(object, format = c("data.table"), ...)

\S4method{getComputeResult}{ComputeResult}(object, format = c("data.table", "igraph"), ...)

\S4method{getComputeResult}{CorrelationResult}(
  object,
  format = c("data.table", "igraph"),
  correlationCoefThreshold = 0.5,
  pValueThreshold = 0.05
)

\S4method{getComputeResult}{DifferentialAbundanceResult}(object, format = c("data.table"))
}
\arguments{
\item{object}{A Microbiome Dataset}

\item{format}{The format of the compute result. Currently only "data.table" and "igraph" are supported.}

\item{...}{additional arguments passed to getComputeResult method of the subclasses of ComputeResult}

\item{correlationCoefThreshold}{threshold to filter edges by correlation coefficient.
Edges with correlation coefficients below this threshold will be removed. Default is .5}

\item{pValueThreshold}{threshold to filter edges by p-value. Edges with p-values above this threshold will be removed. Default is .05}
}
\value{
The compute result in the specified format
}
\description{
Get the compute result from a Microbiome Dataset in a particular format.
Some formats may not be supported for all compute results.
}
