setClassUnion("missingOrNULL", c("missing", "NULL"))

## Helper functions
## will relabel values in X to groupA and groupB based on predicate 
## functions groupAPredicate and groupBPredicate. If groupBPredicate is
## not provided, groupB will be the complement of groupA.
assignToBinaryGroups <- function(x, groupAPredicate, groupBPredicate = NULL) {
    if (!inherits(groupAPredicate, "function")) {
        stop("Argument 'groupAPredicate' must be a function")
    }
    if (!is.null(groupBPredicate)) {
        if (!inherits(groupBPredicate, "function")) {
            stop("Argument 'groupBPredicate' must be a function")
        }
    }
    
    groupAIndexes <- which(groupAPredicate(x))
    if (!is.null(groupBPredicate)) {
        groupBIndexes <- which(groupBPredicate(x))
    } else {
        groupBIndexes <- which(!groupAPredicate(x))
    }

    if (length(intersect(groupAIndexes, groupBIndexes)) > 0) {
        stop("Arguments 'groupA' and 'groupB' must be functions which do not both return true for any given value. They must be made mutually exclusive.")
    }

    x[groupAIndexes] <- "groupA"
    x[groupBIndexes] <- "groupB"

    return(x)
}

#' @importFrom veupathUtils VariableSpec
buildBinaryComparator <- function(covariate, groupAValue, groupBValue) {

    binA <- veupathUtils::Bin(binLabel=groupAValue)
    binB <- veupathUtils::Bin(binLabel=groupBValue)

    groupABins <- veupathUtils::BinList(S4Vectors::SimpleList(c(binA)))
    groupBBins <- veupathUtils::BinList(S4Vectors::SimpleList(c(binB)))

    comparatorVariable <- microbiomeComputations::Comparator(
                            variable = veupathUtils::VariableMetadata(
                                variableSpec = VariableSpec(
                                    variableId = covariate,
                                    entityId = ''
                                ),
                                dataShape = veupathUtils::DataShape(value="CATEGORICAL")
                            ),
                            groupA = groupABins,
                            groupB = groupBBins
    )

    return(comparatorVariable)
}

#### NOTE: were leaving it to microbiomeComputations to dispatch to DESeq2 or Maaslin2 methods 
#### and validate the object we pass to them makes sense. Its possible to pass an AbundanceCollection here
#### w method 'DESeq2'... but microbiomeComputations doesn't currently support that and will complain.
#### I think thats ok, and better than duplicating or complicating logic here.

#' Differential abundance
#'
#' This function returns the fold change and associated p value for
#' a differential abundance analysis. It is useful for finding taxa with
#' an abundance that strongly differs between two groups of samples,
#' such as finding taxa with abundance that differs between skin and saliva samples.
#' This function allows one to recreate the results of differential abundance analysis
#' from MicrobiomeDB.org. However, we recognize that this function makes some assumptions 
#' about the data which may not be valid in other contexts. For better support of 
#' longitudinal studies or metabolomic data, for example, please see our wrapper/ helper methods 
#' for Maaslin2 (\code{MicrobiomeDB::Maaslin2}) and DESeq2 (\code{DESeqDataSetFromCollection}).
#' 
#' @param data AbundanceData object
#' @param covariate character vector giving the name of a metadata variable of interest. If this 
#' variable has only two values, you do not need to provide functions for arguments `groupA` and `groupB`.
#' @param groupA A function that takes a vector of values and returns TRUE or FALSE for each value. This will 
#' be used to assign samples to groupA.
#' @param groupB A function that takes a vector of values and returns TRUE or FALSE for each value. This will 
#' be used to assign samples to groupB. If not provided, groupB will be the complement of groupA.
#' @param method string defining the the differential abundance method. Accepted values are 'DESeq2' and 'Maaslin2'. 
#' Default is 'Maaslin2', as 'DESeq2' only supports counts.
#' @param verbose boolean indicating if timed logging is desired
#' @return ComputeResult object
#' @rdname differentialAbundance-methods
#' @importFrom microbiomeComputations differentialAbundance Comparator
#' @export
setGeneric("differentialAbundance", 
function(data, covariate, groupA, groupB, method = c("Maaslin2", "DESeq2"), verbose = c(TRUE, FALSE)) {
    standardGeneric("differentialAbundance")
}, signature = c("data", "covariate", "groupA", "groupB"))

#' @rdname differentialAbundance-methods
#' @aliases differentialAbundance,CollectionWithMetadata,character,missingOrNULL,missingOrNULL-method
setMethod("differentialAbundance", signature("CollectionWithMetadata", "character", "missingOrNULL", "missingOrNULL"), 
function(data, covariate, groupA, groupB, method = c("Maaslin2", "DESeq2"), verbose = c(TRUE, FALSE)) {
    verbose <- veupathUtils::matchArg(verbose)

    if (data.table::uniqueN(data@sampleMetadata@data[[covariate]]) < 2) {
        stop("Argument 'covariate' must have at least two unique values")
    } else if (data.table::uniqueN(data@sampleMetadata@data[[covariate]]) > 2) {
        stop("Argument 'covariate' must have exactly two unique values if no 'groupA' is provided.")
    }

    groupALabel <- unique(data@sampleMetadata@data[[covariate]])[1]
    groupBLabel <- unique(data@sampleMetadata@data[[covariate]])[2]
    comparator <- buildBinaryComparator(covariate, groupALabel, groupBLabel)
    
    return(microbiomeComputations::differentialAbundance(data, comparator = comparator, method = method, verbose = verbose))
})

#' @rdname differentialAbundance-methods
#' @aliases differentialAbundance,CollectionWithMetadata,character,function,missingOrNULL-method
setMethod("differentialAbundance", signature("CollectionWithMetadata", "character", "function", "missingOrNULL"),
function(data, covariate, groupA, groupB, method = c("Maaslin2", "DESeq2"), verbose = c(TRUE, FALSE)) {
    verbose <- veupathUtils::matchArg(verbose)

    if (data.table::uniqueN(data@sampleMetadata@data[[covariate]]) < 2) {
        stop("Argument 'covariate' must have at least two unique values")
    }

    data@sampleMetadata@data[[covariate]] <- assignToBinaryGroups(data@sampleMetadata@data[[covariate]], groupA, NULL)
    comparator <- buildBinaryComparator(covariate, 'groupA', 'groupB')

    return(microbiomeComputations::differentialAbundance(data, comparator = comparator, method = method, verbose = verbose))
})

#' @rdname differentialAbundance-methods
#' @aliases differentialAbundance,CollectionWithMetadata,character,function,function-method
setMethod("differentialAbundance", signature("CollectionWithMetadata", "character", "function", "function"),
function(data, covariate, groupA, groupB, method = c("Maaslin2", "DESeq2"), verbose = c(TRUE, FALSE)) {
    verbose <- veupathUtils::matchArg(verbose)
    
    if (data.table::uniqueN(data@sampleMetadata@data[[covariate]]) < 2) {
        stop("Argument 'covariate' must have at least two unique values")
    }

    data@sampleMetadata@data[[covariate]] <- assignToBinaryGroups(data@sampleMetadata@data[[covariate]], groupA, groupB)
    comparator <- buildBinaryComparator(covariate, 'groupA', 'groupB')

    ## microbiomeComputations will remove for us rows not in either group, and provide validation
    return(microbiomeComputations::differentialAbundance(data, comparator = comparator, method = method, verbose = verbose))
})

#### NOTE: While i think its important for people to be able to recreate the computes from the site, and tried to make
#### that easier w the differentialAbundance wrapper method above, I also recognize that our
#### microbiomeComputations implementation is limited. So these are methods that should let people pass our
#### objects more-or-less directly to Maaslin2 and DESeq2, using all the same arguments they would usually be able to.

#' Maaslin2 from a Collection
#' 
#' MaAsLin2 finds associations between microbiome 
#' meta-omics features and complex metadata in population-scale 
#' epidemiological studies. The software includes multiple 
#' analysis methods (including support for multiple covariates and repeated measures)
#' filtering, normalization, and transform options to customize analysis for your specific study.
#' 
#' @param data a CollectionWithMetadata
#' @param verbose boolean indicating if timed logging is desired
#' @param ... additional arguments to pass to Maaslin2::Maaslin2
#' @importFrom Maaslin2 Maaslin2
#' @import data.table
#' @rdname Maaslin2
#' @export
setGeneric("Maaslin2", function(data, verbose = c(TRUE,FALSE), ...) standardGeneric("Maaslin2"), signature = c("data"))

#' @rdname Maaslin2
#' @aliases Maaslin2,CollectionWithMetadata-method
setMethod("Maaslin2", signature("CollectionWithMetadata"), function(data, verbose = c(TRUE,FALSE), ...) {
    verbose <- veupathUtils::matchArg(verbose)
    
    recordIdColumn <- data@recordIdColumn
    ancestorIdColumns <- data@ancestorIdColumns
    allIdColumns <- c(recordIdColumn, ancestorIdColumns)
    sampleMetadata <- veupathUtils::getSampleMetadata(data)
    abundances <- microbiomeComputations::getAbundances(data)
    print(class(abundances))

    # remove id columns and any columns that are all 0s.
    cleanedData <- purrr::discard(abundances[, -allIdColumns, with=FALSE], function(col) {identical(union(unique(col), c(0, NA)), c(0, NA))})
    rownames(cleanedData) <- abundances[[recordIdColumn]]
    rownames(sampleMetadata) <- sampleMetadata[[recordIdColumn]]

    maaslinOutput <- Maaslin2::Maaslin2(
        input_data = cleanedData, 
        input_metadata = sampleMetadata,
        ...)

    return(maaslinOutput)
})

#' DESeqDataSet object from a Collection
#'
#' \code{DESeqDataSet} is a subclass of \code{RangedSummarizedExperiment},
#' used to store the input values, intermediate calculations and results of an
#' analysis of differential expression.  The \code{DESeqDataSet} class
#' enforces non-negative integer values in the "counts" matrix stored as
#' the first element in the assay list.
#' In addition, a formula which specifies the design of the experiment must be provided.
#' 
#' @param data AbsoluteAbundanceData object
#' @param verbose boolean indicating if timed logging is desired
#' @param ... additional arguments passed to DESeq2::DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom microbiomeComputations AbsoluteAbundanceData
#' @rdname DESeqDataSetFromCollection
#' @export
setGeneric("DESeqDataSetFromCollection", function(data, verbose = c(TRUE,FALSE), ...) {
    standardGeneric("DESeqDataSetFromCollection")
}, signature = c("data"))

#' @rdname DESeqDataSetFromCollection
#' @aliases DESeqDataSetFromCollection,AbsoluteAbundanceData-method
setMethod("DESeqDataSetFromCollection", signature("AbsoluteAbundanceData"), function(data, verbose = c(TRUE,FALSE), ...) {
    verbose <- veupathUtils::matchArg(verbose)
    
    recordIdColumn <- data@recordIdColumn
    ancestorIdColumns <- data@ancestorIdColumns
    allIdColumns <- c(recordIdColumn, ancestorIdColumns)
    sampleMetadata <- veupathUtils::getSampleMetadata(data)
    abundances <- microbiomeComputations::getAbundances(data, verbose = verbose)

    # First, remove id columns and any columns that are all 0s.
    cleanedData <- purrr::discard(abundances[, -..allIdColumns], function(col) {identical(union(unique(col), c(0, NA)), c(0, NA))})
    # Next, transpose abundance data to get a counts matrix with taxa as rows and samples as columns
    counts <- data.table::transpose(cleanedData)
    rownames(counts) <- names(cleanedData)
    colnames(counts) <- abundances[[recordIdColumn]]

    # Then, format metadata. Recall samples are rows and variables are columns
    rownames(sampleMetadata) <- sampleMetadata[[recordIdColumn]]

    # Finally, check to ensure samples are in the same order in counts and metadata. DESeq
    # expects the order to match, and will not perform this check.
    if (!identical(rownames(sampleMetadata), colnames(counts))){
        # Reorder sampleMetadata to match counts
        veupathUtils::logWithTime("Sample order differs between data and metadata. Reordering data based on the metadata sample order.", verbose)
        data.table::setcolorder(counts, rownames(sampleMetadata))
    }

    # Create DESeqDataSet (dds)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                            colData = sampleMetadata,
                                            ...)

    return(dds)
})