#' Get Microbiome Dataset Compute Result
#' 
#' Get the compute result from a Microbiome Dataset in a particular format.
#' Some formats may not be supported for all compute results.
#' 
#' @examples 
#' correlationOutput <- MicrobiomeDB::selfCorrelation(
#'      getCollection(
#'		microbiomeData::DiabImmune, 
#'		"16S (V4) Genus (Relative taxonomic abundance analysis)"), 
#'      method='spearman', 
#'      verbose=FALSE
#' )
#' correlationDT <- getComputeResult(correlationOutput, "data.table")
#' correlationIGraph <- getComputeResult(correlationOutput, "igraph")
#' @param object A Microbiome Dataset
#' @param format The format of the compute result. Currently only "data.table" and "igraph" are supported.
#' @param ... additional arguments passed to getComputeResult method of the subclasses of ComputeResult
#' @return The compute result in the specified format
#' @importFrom veupathUtils matchArg ComputeResult
#' @export
#' @rdname getComputeResult
setGeneric("getComputeResult", function(object, format = c("data.table"), ...) standardGeneric("getComputeResult"))


#' @rdname getComputeResult
#' @aliases getComputeResult,ComputeResult-method
setMethod("getComputeResult", "ComputeResult", function(object, format = c("data.table", "igraph"), ...) {
    format <- veupathUtils::matchArg(format)

    if (!!length(object@statistics)) {
        return(getComputeResult(object@statistics, format, ...))
    } else {
        if (format == "igraph") {
            stop("igraph not yet supported")
        }
    }

    dt <- data.table::setDT(object@data)

    return(dt)  
})

#' @importFrom veupathUtils CorrelationResult
#' @rdname getComputeResult
#' @param correlationCoefThreshold threshold to filter edges by correlation coefficient. 
#' Edges with correlation coefficients below this threshold will be removed. Default is .5
#' @param pValueThreshold threshold to filter edges by p-value. Edges with p-values above this threshold will be removed. Default is .05
#' @aliases getComputeResult,CorrelationResult-method
#' @importFrom igraph graph_from_data_frame
setMethod("getComputeResult", "CorrelationResult", function(object, format = c("data.table", "igraph"), correlationCoefThreshold = .5, pValueThreshold = .05) {
    format <- veupathUtils::matchArg(format)

    result <- data.table::setDT(object@statistics)

    result <- result[result$correlationCoef >= correlationCoefThreshold & result$pValue <= pValueThreshold, ]

    if (format == "igraph") {
        result <- igraph::graph_from_data_frame(result)
    }

    return(result)  
})

#' @importFrom microbiomeComputations DifferentialAbundanceResult
#' @rdname getComputeResult
#' @aliases getComputeResult,DifferentialAbundanceResult-method
setMethod("getComputeResult", "DifferentialAbundanceResult", function(object, format = c("data.table")) {
    format <- veupathUtils::matchArg(format) 
    return(data.table::setDT(object@statistics))
})

#' @importFrom veupathUtils getSampleMetadata
#' @importFrom veupathUtils getSampleMetadataIdColumns
mergeComputeResultAndMetadata <- function(computeResult, dataset, metadataVariables) {
    dt <- getComputeResult(computeResult, "data.table")
    metadata <- veupathUtils::getSampleMetadata(dataset, includeIds = TRUE, metadataVariables = metadataVariables)

    metadataIdColumns <- veupathUtils::getSampleMetadataIdColumns(dataset)
    dt <- merge(dt, metadata, by = metadataIdColumns, all.x = TRUE)

    return(dt)
}

#' Get Microbiome Dataset Compute Result With Metadata
#' 
#' Get the compute result from a Microbiome Dataset in a particular format with metadata.
#' 
#' @examples 
#' alphaDivOutput <- MicrobiomeDB::alphaDiv(
#'      getCollection(
#'		microbiomeData::DiabImmune, 
#'		"16S (V4) Genus (Relative taxonomic abundance analysis)"), 
#'      method='shannon', 
#'      verbose=FALSE
#' )
#' 
#' alphaDivDT <- getComputeResultWithMetadata(
#'      alphaDivOutput, 
#'      microbiomeData::DiabImmune, 
#'      metadataVariables = c('country', 'delivery_mode')
#' )
#' @param object A Microbiome Dataset
#' @param dataset The MbioDataset, AbundanceData or Collection object from which the compute result was obtained.
#' @param format The format you want the compute result in. Currently only "data.table" is supported.
#' @param metadataVariables The metadata variables to include in the compute result. If NULL, no metadata variables will be included.
#' @return The compute result in the specified format
#' @export
#' @rdname getComputeResultWithMetadata
setGeneric("getComputeResultWithMetadata", 
function(object, dataset, format = c("data.table"), metadataVariables = NULL) 
    standardGeneric("getComputeResultWithMetadata"), 
    signature = c("object", "dataset")
)

#' @rdname getComputeResultWithMetadata
#' @aliases getComputeResultWithMetadata,ComputeResult,MbioDataset-method
setMethod("getComputeResultWithMetadata", signature = c("ComputeResult", "MbioDataset"), 
function(object, dataset = NULL, format = c("data.table"), metadataVariables = NULL) {
    format <- veupathUtils::matchArg(format)
    dt <- mergeComputeResultAndMetadata(object, dataset, metadataVariables)

    return(dt)
})

#' @rdname getComputeResultWithMetadata
#' @aliases getComputeResultWithMetadata,ComputeResult,Collection-method
setMethod("getComputeResultWithMetadata", signature = c("ComputeResult", "Collection"), 
function(object, dataset = NULL, format = c("data.table"), metadataVariables = NULL) {
    format <- veupathUtils::matchArg(format)
    dt <- mergeComputeResultAndMetadata(object, dataset, metadataVariables)

    return(dt)
})

#' @rdname getComputeResultWithMetadata
#' @aliases getComputeResultWithMetadata,ComputeResult,AbundanceData-method
setMethod("getComputeResultWithMetadata", signature = c("ComputeResult", "AbundanceData"), 
function(object, dataset = NULL, format = c("data.table"), metadataVariables = NULL) {
    format <- veupathUtils::matchArg(format)
    dt <- mergeComputeResultAndMetadata(object, dataset, metadataVariables)

    return(dt)
})

#' Correlation Network Visualization
#' 
#' Visualize a correlation result as a network
#' 
#' @examples 
#' correlationOutput <- MicrobiomeDB::correlation(
#'      getCollection(
#'		microbiomeData::DiabImmune, 
#'		"16S (V4) Genus (Relative taxonomic abundance analysis)", 
#'		continuousMetadataOnly = TRUE), 
#'      method='spearman', 
#'      verbose=FALSE
#' )
#' correlationNetwork(correlationOutput) ## renders html widget
#' @param object A ComputeResult or data.frame
#' @param correlationCoefThreshold threshold to filter edges by correlation coefficient.
#' Edges with correlation coefficients below this threshold will be removed. Default is .5
#' @param pValueThreshold threshold to filter edges by p-value. Edges with p-values above this threshold will be removed. Default is .05
#' @param ... additional arguments
#' @export
#' @rdname correlationNetwork
setGeneric("correlationNetwork", function(object, correlationCoefThreshold = .5, pValueThreshold = .05, ...) standardGeneric("correlationNetwork"))

#' @rdname correlationNetwork
#' @aliases correlationNetwork,ComputeResult-method
setMethod("correlationNetwork", "ComputeResult", function(object, correlationCoefThreshold = .5, pValueThreshold = .05, ...) {
    if (!length(object@statistics)) {
        stop("ComputeResult has no statistics")
    }
    if (!inherits(object@statistics, "CorrelationResult")) {
        stop("ComputeResult statistics must be a CorrelationResult")
    }

    edgeList <- getComputeResult(object, correlationCoefThreshold = correlationCoefThreshold, pValueThreshold = pValueThreshold)
    isBipartite <- FALSE
    if (object@computationDetails != "selfCorrelation") {
        isBipartite <- TRUE
    }

    return(suppressWarnings(correlationNetwork(edgeList, bipartiteNetwork = isBipartite)))
})

#' @rdname correlationNetwork
#' @param bipartiteNetwork Should the network use a bipartite or unipartite layout? Defaults to unipartite.
#' @aliases correlationNetwork,data.frame-method
#' @importFrom corGraph bipartiteNetwork
#' @importFrom corGraph unipartiteNetwork
setMethod("correlationNetwork", "data.frame", function(
    object, 
    correlationCoefThreshold = .5, 
    pValueThreshold = .05, 
    bipartiteNetwork = c(FALSE, TRUE)
) {
    bipartiteNetwork <- veupathUtils::matchArg(bipartiteNetwork)

    warning("data.frame input assumes the columns are in the following order: source, target, correlationCoef, pValue.")
    names(object) <- c("source", "target", "value", "pValue")

    object <- object[object$value >= correlationCoefThreshold & object$pValue <= pValueThreshold, ]

    sources <- unique(object$source)
    targets <- unique(object$target)

    #are all sources different from targets?
    if (bipartiteNetwork) {
        net <- corGraph::bipartiteNetwork(object)
    } else {
        net <- corGraph::unipartiteNetwork(object)
    }

    return(net)
})
