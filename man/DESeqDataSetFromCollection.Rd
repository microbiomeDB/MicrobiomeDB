% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-differentialAbundance.R
\name{DESeqDataSetFromCollection}
\alias{DESeqDataSetFromCollection}
\alias{DESeqDataSetFromCollection,AbsoluteAbundanceData-method}
\title{DESeqDataSet object from a Collection}
\usage{
DESeqDataSetFromCollection(data, verbose = c(TRUE, FALSE), ...)

\S4method{DESeqDataSetFromCollection}{AbsoluteAbundanceData}(data, verbose = c(TRUE, FALSE), ...)
}
\arguments{
\item{data}{AbsoluteAbundanceData object}

\item{verbose}{boolean indicating if timed logging is desired}

\item{...}{additional arguments passed to DESeq2::DESeqDataSetFromMatrix}
}
\description{
\code{DESeqDataSet} is a subclass of \code{RangedSummarizedExperiment},
used to store the input values, intermediate calculations and results of an
analysis of differential expression.  The \code{DESeqDataSet} class
enforces non-negative integer values in the "counts" matrix stored as
the first element in the assay list.
In addition, a formula which specifies the design of the experiment must be provided.
}
