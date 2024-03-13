#' Get Microbiome Dataset Compute Result
#' 
#' Get the compute result from a Microbiome Dataset in a particular format.
#' Some formats may not be supported for all compute results.
#' @param object A Microbiome Dataset
#' @param format The format of the compute result. Currently only "data.table" and "igraph" are supported.
#' @return The compute result in the specified format
#' @importFrom veupathUtils matchArg ComputeResult
#' @export
#' @rdname getComputeResult
setGeneric("getComputeResult", function(object, format = c("data.table"), ...) standardGeneric("getComputeResult"))


#' @rdname getComputeResult
#' @aliases getComputeResult,ComputeResult-method
setMethod("getComputeResult", "ComputeResult", function(object, format = c("data.table", "igraph")) {
    format <- veupathUtils::matchArg(format)

    if (!!length(object@statistics)) {
        return(getComputeResult(object@statistics, format))
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
#' @aliases getComputeResult,CorrelationResult-method
setMethod("getComputeResult", "CorrelationResult", function(object, format = c("data.table", "igraph")) {
    format <- veupathUtils::matchArg(format)

    result <- data.table::setDT(object@statistics)

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

#' @importFrom microbiomeData getSampleMetadata
#' @importFrom microbiomeData getSampleMetadataIdColumns
mergeComputeResultAndMetadata <- function(computeResult, dataset, metadataVariables) {
    dt <- getComputeResult(computeResult, "data.table")
    metadata <- microbiomeData::getSampleMetadata(dataset, includeIds = TRUE, metadataVariables = metadataVariables)

    metadataIdColumns <- microbiomeData::getSampleMetadataIdColumns(dataset)
    print(metadataIdColumns)
    print(head(names(dt)))
    print(head(names(metadata)))
    dt <- merge(dt, metadata, by = metadataIdColumns, all.x = TRUE)

    return(dt)
}

#' Get Microbiome Dataset Compute Result With Metadata
#' 
#' Get the compute result from a Microbiome Dataset in a particular format with metadata.
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