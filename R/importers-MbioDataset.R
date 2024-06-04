
buildCollectionFromTreeSE <- function(
    collectionName = list(assayDataName = NULL, rowDataColumnName = NULL), 
    rowData, 
    assayData,
    normalizationMethod = c("TSS", "none"),
    verbose = c(TRUE, FALSE)
) {
    verbose <- veupathUtils::matchArg(verbose)

    assayDataName <- collectionName$assayDataName
    rowDataColumnName <- collectionName$rowDataColumnName

    if (is.null(assayDataName) || is.null(rowDataColumnName)) {
        stop("Must specify both assayDataName and rowDataColumnName as named elements of the collectionName list argument")
    }

    assayDT <- as.data.frame.matrix(assayData, col.names = colnames(assayData), row.names = row.names(assayData))
    dt <- data.table::as.data.table(merge(assayDT, rowData[rowDataColumnName], by = 0))
    dt$Row.names <- NULL

    recordIDs <- names(dt)[names(dt) != rowDataColumnName]
    dt <- dt[, lapply(.SD, sum, na.rm=TRUE), by=rowDataColumnName]
    dt <- data.table::transpose(dt, make.names=rowDataColumnName)

    # if this does grow into other methods, the normalization step could be factored out probably
    if (normalizationMethod == "TSS") {
        dt <- dt / rowSums(dt)
    }

    dt$recordIDs <- recordIDs

    recordIdColumn <- 'recordIDs'
    ancestorIdColumns <- NULL
    collectionName <- paste0(assayDataName, ": ", rowDataColumnName)
    if (normalizationMethod != "none") {
        collectionName <- paste0(collectionName, " (", normalizationMethod, " normalized)")
    }

    collection <- Collection(
        data = dt,
        recordIdColumn = recordIdColumn,
        ancestorIdColumns = ancestorIdColumns,
        collectionName = collectionName
    )

    return(collection)
}

#' Import TreeSummarizedExperiment
#' 
#' Import data from TreeSummarizedExperiment to MbioDataset.
#' There is some loss of granularity in this process. It results
#' in a simpler and more performant object which is compliant
#' with the MicrobiomeDB infrastructure.
#' 
#' @param data A TreeSummarizedExperiment
#' @param normalizationMethod Normalization method to use on they assay data. Options are "none" and "TSS".
#' Applying TSS normalization to absolute taxonomic abundances produces relative taxonomic abundances. Default is "TSS".
#' @param keepRawValues Keep the raw assay values as well as the normalized values.
#' @param verbose Print messages
#' @return A MbioDataset
#' @export
importTreeSummarizedExperiment <- function(data, normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE)) {

    normalizationMethod <- veupathUtils::matchArg(normalizationMethod)
    keepRawValues <- veupathUtils::matchArg(keepRawValues)
    verbose <- veupathUtils::matchArg(verbose)

    if (keepRawValues == FALSE && normalizationMethod == "none") {
        stop("keepRawValues must be TRUE when normalizationMethod is 'none'")
    }

    # figure out what all assays we have, and what data is available per assay
    # these will become collections in the MbioDataset
    collectionsDTList <- lapply(names(data@assays), function(x) {
        data.table::data.table(assayDataName = x, rowDataColumnName = row.names(data@assays[[x]]))
    })
    collectionsDT <- purrr::reduce(collectionsDT, rbind)

    if (keepRawValues) {
        # call buildCollectionFromTreeSE for each column of each assay data, w normalization 'none'
        collectionsByAssayList <- apply(collectionsDT, 1, function(x) {
            buildCollectionFromTreeSE(
                collectionName = as.list(x, keep.names=TRUE), 
                rowData = data@rowData, 
                assayData = data@assays[[x$assayDataName]], 
                normalizationMethod = "none",
                verbose = verbose
            )
        })
        rawCollectionsList <- purrr::reduce(collectionsByAssayList, c)
    }

    if (normalizationMethod != "none") {
        collectionsByAssayList <- apply(collectionsDT, 1, function(x) {
            buildCollectionFromTreeSE(
                collectionName = as.list(x, keep.names=TRUE), 
                rowData = data@rowData, 
                assayData = data@assays[[x$assayDataName]], 
                normalizationMethod = normalizationMethod, 
                verbose = verbose
            )
        })
        normalizedCollectionsList <- purrr::reduce(collectionsByAssayList, c)
    }
    
    collectionsList <- c(rawCollectionsList, normalizedCollectionsList)

    ## TODO add some sort of check and set this conditionally
    imputeZero <- FALSE

    # build and validate MbioDataset, colData becomes sampleMetadata
    mbioDataset <- MbioDataset(
        collections = collectionsList, 
        metadata = SampleMetadata(data@colData),
        imputeZero = imputeZero)

    # return a MbioDataset
    return(mbioDataset)
}

## lean on miaverse to import biom, phyloseq, csv, etc
## TODO do these also needs args about relative abundances? id think so..

#' Import HUMAnN data
#' 
#' Import data from HUMAnN results to MbioDataset. There is
#' some loss of granularity in this process. It results
#' in a simpler and more performant object which is compliant
#' with the MicrobiomeDB infrastructure. See \code{mia::importHUMAnN}
#' for documentation.
#' 
#' @param normalizationMethod Normalization method to use on they assay data. Options are "none" and "TSS".
#' Applying TSS normalization to absolute taxonomic abundances produces relative taxonomic abundances. Default is "TSS".
#' @param keepRawValues Keep the raw assay values as well as the normalized values.
#' @param verbose Print messages
#' @param ... Arguments to pass to mia::importHUMAnN
#' @return A MbioDataset
#' @export
#' @importFrom mia importHUMAnN
importHUMAnN <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    treeSE <- mia::importHUMAnN(...)

    mbioDataset <- importTreeSummarizedExperiment(treeSE, normalizationMethod = normalizationMethod, keepRawValues = keepRawValues, verbose = verbose)
    return(mbioDataset)
}

#' Import MetaPhlAn data
#' 
#' Import data from MetaPhlAn results to MbioDataset. There is
#' some loss of granularity in this process. It results
#' in a simpler and more performant object which is compliant
#' with the MicrobiomeDB infrastructure. See \code{mia::importMetaPhlAn}
#' for documentation.
#' 
#' @param normalizationMethod Normalization method to use on they assay data. Options are "none" and "TSS".
#' Applying TSS normalization to absolute taxonomic abundances produces relative taxonomic abundances. Default is "TSS".
#' @param keepRawValues Keep the raw assay values as well as the normalized values.
#' @param verbose Print messages
#' @param ... Arguments to pass to mia::importMetaPhlAn
#' @return A MbioDataset
#' @export
#' @importFrom mia importMetaPhlAn
importMetaPhlAn <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    treeSE <- mia::importMetaPhlAn(...)

    mbioDataset <- importTreeSummarizedExperiment(treeSE, normalizationMethod = normalizationMethod, keepRawValues = keepRawValues, verbose = verbose)
    return(mbioDataset)
}

#' Import MOTHUR data
#' 
#' Import data from MOTHUR results to MbioDataset. There is
#' some loss of granularity in this process. It results
#' in a simpler and more performant object which is compliant
#' with the MicrobiomeDB infrastructure. See \code{mia::importMothur}
#' for documentation.
#' 
#' @param normalizationMethod Normalization method to use on they assay data. Options are "none" and "TSS".
#' Applying TSS normalization to absolute taxonomic abundances produces relative taxonomic abundances. Default is "TSS".
#' @param keepRawValues Keep the raw assay values as well as the normalized values.
#' @param verbose Print messages
#' @param ... Arguments to pass to mia::importMothur
#' @return A MbioDataset
#' @export
#' @importFrom mia importMothur
importMothur <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    treeSE <- mia::importMothur(...)

    mbioDataset <- importTreeSummarizedExperiment(treeSE, normalizationMethod = normalizationMethod, keepRawValues = keepRawValues, verbose = verbose)
    return(mbioDataset)
}

#' Import QIIME2 data
#' 
#' Import data from QIIME2 results to MbioDataset. There is
#' some loss of granularity in this process. It results
#' in a simpler and more performant object which is compliant
#' with the MicrobiomeDB infrastructure. See \code{mia::importQIIME2}
#' for documentation.
#' 
#' @param normalizationMethod Normalization method to use on they assay data. Options are "none" and "TSS".
#' Applying TSS normalization to absolute taxonomic abundances produces relative taxonomic abundances. Default is "TSS".
#' @param keepRawValues Keep the raw assay values as well as the normalized values.
#' @param verbose Print messages
#' @param ... Arguments to pass to mia::importQIIME2
#' @return A MbioDataset
#' @export
#' @importFrom mia importQIIME2
importQIIME2 <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    treeSE <- mia::importQIIME2(...)

    mbioDataset <- importTreeSummarizedExperiment(treeSE, normalizationMethod = normalizationMethod, keepRawValues = keepRawValues, verbose = verbose)
    return(mbioDataset)
}

#' Import BIOM data
#' 
#' Import data from BIOM results to MbioDataset. There is
#' some loss of granularity in this process. It results
#' in a simpler and more performant object which is compliant
#' with the MicrobiomeDB infrastructure. See \code{mia::makeTreeSEFromBiom}
#' for documentation.
#' 
#' @param normalizationMethod Normalization method to use on they assay data. Options are "none" and "TSS".
#' Applying TSS normalization to absolute taxonomic abundances produces relative taxonomic abundances. Default is "TSS".
#' @param keepRawValues Keep the raw assay values as well as the normalized values.
#' @param verbose Print messages
#' @param ... Arguments to pass to mia::makeTreeSEFromBiom
#' @return A MbioDataset
#' @export
#' @importFrom mia makeTreeSEFromBiom
importBIOM <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    treeSE <- mia::makeTreeSEFromBiom(...)

    mbioDataset <- importTreeSummarizedExperiment(treeSE, normalizationMethod = normalizationMethod, keepRawValues = keepRawValues, verbose = verbose)
    return(mbioDataset)
}

#' Import DADA2 data
#' 
#' Import data from DADA2 results to MbioDataset. There is
#' some loss of granularity in this process. It results
#' in a simpler and more performant object which is compliant
#' with the MicrobiomeDB infrastructure. See \code{mia::makeTreeSEFromDADA2}
#' for documentation.
#' 
#' @param normalizationMethod Normalization method to use on they assay data. Options are "none" and "TSS".
#' Applying TSS normalization to absolute taxonomic abundances produces relative taxonomic abundances. Default is "TSS".
#' @param keepRawValues Keep the raw assay values as well as the normalized values.
#' @param verbose Print messages
#' @param ... Arguments to pass to mia::makeTreeSEFromDADA2
#' @return A MbioDataset
#' @export
#' @importFrom mia makeTreeSEFromDADA2
importDADA2 <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    treeSE <- mia::makeTreeSEFromDADA2(...)

    mbioDataset <- importTreeSummarizedExperiment(treeSE, normalizationMethod = normalizationMethod, keepRawValues = keepRawValues, verbose = verbose)
    return(mbioDataset)
}

#' Import Phyloseq data
#' 
#' Import data from Phyloseq results to MbioDataset. There is
#' some loss of granularity in this process. It results
#' in a simpler and more performant object which is compliant
#' with the MicrobiomeDB infrastructure. See \code{mia::makeTreeSEFromPhyloseq}
#' for documentation.
#' 
#' @param normalizationMethod Normalization method to use on they assay data. Options are "none" and "TSS".
#' Applying TSS normalization to absolute taxonomic abundances produces relative taxonomic abundances. Default is "TSS".
#' @param keepRawValues Keep the raw assay values as well as the normalized values.
#' @param verbose Print messages
#' @param ... Arguments to pass to mia::makeTreeSEFromPhyloseq
#' @return A MbioDataset
#' @export
#' @importFrom mia makeTreeSEFromPhyloseq
importPhyloseq <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    treeSE <- mia::makeTreeSEFromPhyloseq(...)

    mbioDataset <- importTreeSummarizedExperiment(treeSE, normalizationMethod = normalizationMethod, keepRawValues = keepRawValues, verbose = verbose)
    return(mbioDataset)
}