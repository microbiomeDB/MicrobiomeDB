#' Ranked abundance
#'
#' This function returns abundances, ranked by a selected ranking function
#' 
#' @examples
#' rankedAbundOutput <- rankedAbundance(getCollection(microbiomeData::DiabImmune, "16S (V4) Genus"), method = "median")
#' @param data AbundanceData object
#' @param method string defining the ranking strategy by which to order the taxa. Accepted values are 'median','max','q3',and 'variance'. Note that taxa that return a value of 0 for a given method will not be included in the results.
#' @param cutoff integer indicating the maximium number of taxa to be kept after ranking.
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @importFrom microbiomeComputations rankedAbundance
#' @export
#' @rdname rankedAbundance-methods
microbiomeComputations::rankedAbundance

#' Alpha diversity
#'
#' This function returns alpha diversity values for each sample.
#' 
#' @examples
#' alphaDivOutput <- alphaDiv(getCollection(microbiomeData::DiabImmune, "16S (V4) Genus"), method = "shannon")
#' @param data AbundanceData object
#' @param method string defining the the alpha diversity method. Accepted values are 'shannon','simpson', and 'evenness'
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @importFrom microbiomeComputations alphaDiv
#' @export
#' @rdname alphaDiv-methods
microbiomeComputations::alphaDiv

#' Beta diversity
#'
#' This function returns pcoa coordinates calculated from the beta diversity dissimilarity matrix.
#' 
#' @examples 
#' betaDivOutput <- betaDiv(getCollection(microbiomeData::DiabImmune, "16S (V4) Genus"), method = "bray", k = 2)
#' @param data AbundanceData object
#' @param method string defining the the beta diversity dissimilarity method. Accepted values are 'bray','jaccard', and 'jsd'
#' @param k integer determining the number of pcoa axes to return
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @importFrom microbiomeComputations betaDiv
#' @export
#' @rdname betaDiv-methods
microbiomeComputations::betaDiv

#' Correlation
#'
#' This function returns correlation coefficients for variables in one dataset against variables in a second dataset
#' 
#' @examples
#' correlationDT <- correlation(getCollection(microbiomeData::DiabImmune, "16S (V4) Genus"), method = 'spearman', format = 'data.table')
#' correlationOutput <- correlation(getCollection(microbiomeData::DiabImmune, "16S (V4) Genus"), method = 'spearman', format = 'ComputeResult')
#' alsoCorrelationDT <- getComputeResult(correlationOutput, "data.table")
#' @param data1 first dataset. A data.table
#' @param data2 second dataset. A data.table
#' @param method string defining the type of correlation to run. 
#' The currently supported values are specific to the class of data1 and data2.
#' @param format string defining the desired format of the result. 
#' The currently supported values are 'data.table' and 'ComputeResult'.
#' @param verbose boolean indicating if timed logging is desired
#' @param ... additional parameters
#' @return data.frame with correlation coefficients or a ComputeResult object
#' @importFrom veupathUtils correlation
#' @export
#' @rdname correlation-methods
veupathUtils::correlation

#' Self Correlation
#'
#' This function returns correlation coefficients for variables in one AbundanceData object against itself. It generally serves as a 
#' convenience wrapper around veupathUtils::correlation, with the exception that it additionally supports sparcc.
#' 
#' @examples
#' correlationDT <- selfCorrelation(getCollection(microbiomeData::DiabImmune, "16S (V4) Genus"), method = 'sparcc', format = 'data.table')
#' correlationOutput <- selfCorrelation(getCollection(microbiomeData::DiabImmune, "16S (V4) Genus"), method = 'sparcc', format = 'ComputeResult')
#' alsoCorrelationDT <- getComputeResult(correlationOutput, "data.table")
#' @param data An AbundanceData object
#' @param method string defining the type of correlation to run. The currently supported values are 'spearman','pearson' and 'sparcc'
#' @param format string defining the desired format of the result. The currently supported values are 'data.table' and 'ComputeResult'.
#' @param verbose boolean indicating if timed logging is desired
#' @param proportionNonZeroThreshold numeric threshold to filter features by proportion of non-zero values across samples
#' @param varianceThreshold numeric threshold to filter features by variance across samples
#' @param stdDevThreshold numeric threshold to filter features by standard deviation across samples
#' @param ... additional parameters
#' @return ComputeResult object
#' @importFrom veupathUtils selfCorrelation
#' @importFrom microbiomeComputations selfCorrelation
#' @export
#' @rdname selfCorrelation-methods
veupathUtils::selfCorrelation
