
buildCollectionFromTreeSE <- function(
    collectionName = list(assayDataName = NULL, rowDataColumnName = NULL), 
    rowData, # this is a data.frame representing the row data/ tree
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
    ancestorIdColumns <- character(0)
    collectionName <- paste0(assayDataName, ": ", rowDataColumnName)
    if (normalizationMethod != "none") {
        collectionName <- paste0(collectionName, " (", normalizationMethod, " normalized)")
    }

    collection <- veupathUtils::Collection(
        data = dt,
        recordIdColumn = recordIdColumn,
        ancestorIdColumns = ancestorIdColumns,
        name = collectionName
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
#' @importFrom purrr reduce
#' @importFrom data.table data.table
#' @importFrom SummarizedExperiment rowData
#' @import TreeSummarizedExperiment
#' @rdname importTreeSummarizedExperiment
#' @examples 
#' data(GlobalPatterns, package="mia")
#' tse <- GlobalPatterns
#' 
#' ## no normalization, with raw values
#' mbioDataset <- importTreeSummarizedExperiment(
#'      tse, 
#'      normalizationMethod = "none", 
#'      keepRawValues = TRUE, 
#'      verbose = TRUE
#' )
#' 
#' ## TSS normalization, drop raw values
#' mbioDataset <- importTreeSummarizedExperiment(
#'      tse, 
#'      normalizationMethod = "TSS", 
#'      keepRawValues = FALSE, 
#'      verbose = TRUE
#' )
#' 
#' ## TSS normalization, keep raw values
#' mbioDataset <- importTreeSummarizedExperiment(
#'      tse, 
#'      normalizationMethod = "TSS", 
#'      keepRawValues = TRUE, 
#'      verbose = TRUE
#' )
#' 
#' @export
importTreeSummarizedExperiment <- function(data, normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE)) {

    normalizationMethod <- veupathUtils::matchArg(normalizationMethod)
    keepRawValues <- veupathUtils::matchArg(keepRawValues)
    verbose <- veupathUtils::matchArg(verbose)

    if (!inherits(data, "SummarizedExperiment")) {
        stop("data must be or extend a SummarizedExperiment")
    }

    if (keepRawValues == FALSE && normalizationMethod == "none") {
        stop("keepRawValues must be TRUE when normalizationMethod is 'none'")
    }

    # figure out what all assays we have, and what data is available per assay
    # these will become collections in the MbioDataset
    # TODO it looks like rowData is expected to be same across all assays,
    # which is odd to me, but probably means we can simplify this logic
    collectionsDTList <- lapply(names(data@assays@data), function(x) {
        data.table::data.table(assayDataName = x, rowDataColumnName = colnames(SummarizedExperiment::rowData(data)))
    })
    collectionsDT <- purrr::reduce(collectionsDTList, rbind)

    if (nrow(collectionsDT) != 0) {
        rawCollectionsList <- list()
        if (keepRawValues) {
            # call buildCollectionFromTreeSE for each column of each assay data, w normalization 'none'
            collectionsByAssayList <- apply(collectionsDT, 1, function(x) {
                collectionName = as.list(x, keep.names=TRUE);

                buildCollectionFromTreeSE(
                    collectionName = collectionName, 
                    rowData = as.data.frame(SummarizedExperiment::rowData(data)), 
                    assayData = data@assays@data[[collectionName$assayDataName]], 
                    normalizationMethod = "none",
                    verbose = verbose
                )
            })
            rawCollectionsList <- purrr::reduce(collectionsByAssayList, c)
        }

        normalizedCollectionsList <- list()
        if (normalizationMethod != "none") {
            collectionsByAssayList <- apply(collectionsDT, 1, function(x) {
                collectionName = as.list(x, keep.names=TRUE);

                buildCollectionFromTreeSE(
                    collectionName = collectionName,
                    rowData = as.data.frame(SummarizedExperiment::rowData(data)), 
                    assayData = data@assays@data[[collectionName$assayDataName]], 
                    normalizationMethod = normalizationMethod, 
                    verbose = verbose
                )
            })
            normalizedCollectionsList <- purrr::reduce(collectionsByAssayList, c)
        }
    
        collectionsList <- c(rawCollectionsList, normalizedCollectionsList)
    } else {
        collectionsList <- veupathUtils::Collections()
    }
    

    # build and validate MbioDataset, colData becomes sampleMetadata
    colData <- SummarizedExperiment::colData(data)
    metadataDT <- data.table::data.table()
    if (!!length(colData)) {
        metadataDT <- data.table::as.data.table(SummarizedExperiment::colData(data))
        metadataDT$recordIDs <- rownames(SummarizedExperiment::colData(data))
    }
    if (length(metadataDT) == 1) metadataDT <- data.table::data.table()

    mbioDataset <- MbioDataset(
        collections = collectionsList, 
        metadata = SampleMetadata(data = metadataDT, recordIdColumn = "recordIDs")
    )

    # return a MbioDataset
    return(mbioDataset)
}

#' @rdname importTreeSummarizedExperiment
importTreeSE <- importTreeSummarizedExperiment

## lean on miaverse to import biom, phyloseq, csv, etc

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
#' @examples 
#' file_path <- system.file("extdata", "humann_output.tsv", package = "mia")
#' mbioDataset <- importHUMAnN(
#'      normalizationMethod = "none", 
#'      keepRawValues = TRUE, 
#'      verbose = TRUE, 
#'      file_path
#' )
#' @export
#' @importFrom mia importHUMAnN
importHUMAnN <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    .require_package("mia")

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
importMetaPhlAn <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    .require_package("mia")

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
#' @examples
#' counts <- system.file("extdata", "mothur_example.shared", package = "mia")
#' taxa <- system.file("extdata", "mothur_example.cons.taxonomy", package = "mia")
#' meta <- system.file("extdata", "mothur_example.design", package = "mia")
#' 
#' mbioDataset <- importMothur(
#'      normalizationMethod = "none", 
#'      keepRawValues = TRUE, 
#'      verbose = TRUE, 
#'      counts, 
#'      taxa, 
#'      meta
#' )
#' @export
importMothur <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    .require_package("mia")

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
#' @examples
#' featureTableFile <- system.file("extdata", "table.qza", package = "mia")
#' taxonomyTableFile <- system.file("extdata", "taxonomy.qza", package = "mia")
#' 
#' mbioDataset <- importQIIME2(
#'      normalizationMethod = "none", 
#'      keepRawValues = TRUE, 
#'      verbose = TRUE, 
#'      featureTableFile, 
#'      taxonomyTableFile
#' )
#' @export
importQIIME2 <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    .require_package("mia")

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
#' @examples
#' rich_dense_file = system.file("extdata", "rich_dense_otu_table.biom",
#'                                package = "biomformat")
#'
#' mbioDataset <- importBIOM(
#'      normalizationMethod = "none", 
#'      keepRawValues = TRUE, 
#'      verbose = TRUE, 
#'      rich_dense_file
#' )
#' 
#' rich_dense_biom = biomformat::read_biom(rich_dense_file)
#' 
#' mbioDataset <- importBIOM(
#'      normalizationMethod = "none", 
#'      keepRawValues = TRUE, 
#'      verbose = TRUE, 
#'      rich_dense_biom
#' )
#' @export
importBIOM <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    .require_package("mia")
    .require_package("biomformat")

    if (!inherits(..., "biom")) {
        biom <- biomformat::read_biom(...)
        treeSE <- mia::makeTreeSEFromBiom(obj=biom)
    } else {
        treeSE <- mia::makeTreeSEFromBiom(...)
    }
    
    

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
#' @examples
#' fnF <- system.file("extdata", "sam1F.fastq.gz", package="dada2")
#' fnR = system.file("extdata", "sam1R.fastq.gz", package="dada2")
#' dadaF <- dada2::dada(fnF, selfConsist=TRUE)
#' dadaR <- dada2::dada(fnR, selfConsist=TRUE)
#' 
#' mbioDataset <- importDADA2(
#'      normalizationMethod = "none", 
#'      keepRawValues = TRUE, 
#'      verbose = TRUE, 
#'      dadaF, 
#'      fnF, 
#'      dadaR, 
#'      fnR
#' )
#' @export
importDADA2 <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    .require_package("mia")
    .require_package("dada2")

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
#' @examples
#' data(GlobalPatterns, package="phyloseq")
#' 
#' mbioDataset <- importPhyloseq(
#'      normalizationMethod = "none", 
#'      keepRawValues = TRUE, 
#'      verbose = TRUE, 
#'      GlobalPatterns
#' )
#' @export
importPhyloseq <- function(normalizationMethod = c("TSS", "none"), keepRawValues = c(TRUE, FALSE), verbose = c(TRUE, FALSE), ...) {
    .require_package("mia")
    .require_package("phyloseq")

    treeSE <- mia::makeTreeSEFromPhyloseq(...)

    mbioDataset <- importTreeSummarizedExperiment(treeSE, normalizationMethod = normalizationMethod, keepRawValues = keepRawValues, verbose = verbose)
    return(mbioDataset)
}