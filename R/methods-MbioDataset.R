collectionNamesGeneric <- getGeneric("getCollectionNames", "veupathUtils")
#' Get Names of Collections
#' 
#' Get the names of the collections in a MbioDataset object
#' 
#' @examples 
#' getCollectionNames(microbiomeData::DiabImmune)
#' @param object An MbioDataset
#' @return A character vector of collection names
#' @export
#' @importFrom veupathUtils getCollectionNames
setMethod(collectionNamesGeneric, "MbioDataset", function(object) return(unname(getCollectionNames(object@collections))))

metadataVarNamesGeneric <- getGeneric("getMetadataVariableNames", "veupathUtils")
#' Get Variable Names of Metadata
#' 
#' Get the names of the metadata variables in an MbioDataset.
#' 
#' @examples 
#' getMetadataVariableNames(microbiomeData::DiabImmune)
#' @param object An MbioDataset
#' @return a character vector of metadata variable names
#' @export
#' @importFrom veupathUtils getMetadataVariableNames
setMethod(metadataVarNamesGeneric, "MbioDataset", function(object) return(names(object@metadata@data)))

metadataVarSummaryGeneric <- getGeneric("getMetadataVariableSummary", "veupathUtils")
#' Get Summary of Metadata Variables
#' 
#' Get a summary of the requested metadata variable in an MbioDataset.
#' 
#' @examples 
#' getMetadataVariableSummary(microbiomeData::DiabImmune, "age_months")
#' getMetadataVariableSummary(microbiomeData::DiabImmune, "sex")
#' getMetadataVariableSummary(microbiomeData::DiabImmune, "country")
#' @param object An MbioDataset
#' @param variable A character vector representing the name of the metadata variable to summarize
#' @return a table summarizing the values of the requested metadata variable
#' @export
#' @importFrom veupathUtils getMetadataVariableSummary
setMethod(metadataVarSummaryGeneric, "MbioDataset", function(object, variable) {
    if (!variable %in% getMetadataVariableNames(object)) {
        stop("Variable ", variable, " not found in sample metadata. Available variables: ", paste(getMetadataVariableNames(object), collapse = ", "))
    }

    varData <- object@metadata@data[[variable]]

    if (class(varData) %in% c('factor','character')) {
        out <- table(varData)
        dimnames(out) <- unname(dimnames(out))
        return(out)
    } else {
        return(summary(varData))
    }
})

sampleMetadataGeneric <- getGeneric("getSampleMetadata", "veupathUtils")
#' Get data.table of sample metadata from MbioDataset
#'
#' Returns a data.table of sample metadata
#' 
#' @examples
#' getSampleMetadata(microbiomeData::DiabImmune)
#' getSampleMetadata(microbiomeData::DiabImmune, metadataVariables = c("age_months", "sex"))
#' @param object MbioDataset
#' @param asCopy boolean indicating whether to return the data as a copy or by reference
#' @param includeIds boolean indicating whether we should include recordIdColumn and ancestorIdColumns
#' @param metadataVariables The metadata variables to include in the sample metadata. If NULL, all metadata variables will be included.
#' @return data.table of sample metadata
#' @export
setMethod(sampleMetadataGeneric, "MbioDataset", function(object, asCopy = c(TRUE, FALSE), includeIds = c(TRUE, FALSE), metadataVariables = NULL) {
    asCopy <- veupathUtils::matchArg(asCopy)
    includeIds <- veupathUtils::matchArg(includeIds)

    if (!length(object@metadata@data)) return(NULL)

    dt <- object@metadata@data
    allIdColumns <- veupathUtils::getSampleMetadataIdColumns(object)

    # Check that incoming dt meets requirements
    if (!inherits(dt, 'data.table')) {
        data.table::setDT(dt)
    }

    if (asCopy) {
        dt <- data.table::copy(dt)
    }

    if (includeIds && !is.null(metadataVariables)) {
        dt <- dt[, c(allIdColumns, metadataVariables), with = FALSE]
    } else if (!includeIds && !is.null(metadataVariables)) {
        dt <- dt[, metadataVariables, with = FALSE]
    } else if (!includeIds && is.null(metadataVariables)) {
        dt <- dt[, -..allIdColumns]
    }

    return(dt)
})

metadataIdColsGeneric <- getGeneric("getSampleMetadataIdColumns", "veupathUtils")
#' Get Sample Metadata Id Column Names
#' 
#' Get the names of the record and ancestor id columns in the sample metadata of an MbioDataset object.
#' @param object MbioDataset
#' @return a character vector of id column names
setMethod(metadataIdColsGeneric, "MbioDataset", function(object) veupathUtils::getIdColumns(object@metadata))


#' Update Microbiome Dataset Collection Name
#' 
#' Update the name of a collection in the Microbiome Dataset.
#' 
#' @examples
#' myCopyOfDiabImmune <- microbiomeData::DiabImmune
#' myCopyOfDiabImmune <- updateCollectionName(
#'	myCopyOfDiabImmune, 
#'	"16S (V4) Genus (Relative taxonomic abundance analysis)", 
#'	"16S Genus"
#' )
#' @param object A Microbiome Dataset
#' @param oldName The name of the collection to update
#' @param newName The new name of the collection
#' @return A Microbiome Dataset with the updated collection name
#' @rdname updateCollectionName
#' @export
setGeneric("updateCollectionName", function(object, oldName, newName) standardGeneric("updateCollectionName"))

#' @rdname updateCollectionName
#' @aliases updateCollectionName,MbioDataset,character,character-method
setMethod("updateCollectionName", "MbioDataset", function(object, oldName, newName) {
    object@collections[[which(getCollectionNames(object) == oldName)]]@name <- newName
    return(object)
})

#' Get Microbiome Dataset Collection
#' 
#' Get a collection from the Microbiome Dataset. The collection will be returned
#' as an AbundanceData, phyloseq, or Collection object.
#' 
#' @examples
#' genus <- getCollection(
#'      microbiomeData::DiabImmune, 
#'      "16S (V4) Genus (Relative taxonomic abundance analysis)"
#' )
#' 
#' genus_phyloseq <- getCollection(
#'      microbiomeData::DiabImmune, 
#'      "16S (V4) Genus (Relative taxonomic abundance analysis)", 
#'      format = "phyloseq"
#' )
#' 
#' ## to pass to correlation method, we want only continuous metadata
#' genus_continuous <- getCollection(
#'      microbiomeData::DiabImmune, 
#'      "16S (V4) Genus (Relative taxonomic abundance analysis)", 
#'      continuousMetadataOnly = TRUE
#' ) 
#' 
#' ## with no metadata
#' genus_collection <- getCollection(
#'      microbiomeData::DiabImmune, 
#'      "16S (V4) Genus (Relative taxonomic abundance analysis)", 
#'      format = "Collection"
#' )
#' @param object A Microbiome Dataset
#' @param collectionName The name of the collection to return
#' @param format The format of the collection to return. Currently supported options are "AbundanceData", "phyloseq" and "Collection".
#' @param continuousMetadataOnly If TRUE, only continuous metadata will be returned. If FALSE, all metadata will be returned.
#' @return An AbundanceData, phyloseq, or Collection object representing the collection and any associated study metadata
#' @importFrom phyloseq phyloseq
#' @importFrom microbiomeComputations AbundanceData
#' @rdname getCollection
#' @export
setGeneric("getCollection", function(object, collectionName, format = c("AbundanceData", "phyloseq", "Collection"), continuousMetadataOnly = c(FALSE, TRUE)) standardGeneric("getCollection"))

#' @rdname getCollection
#' @aliases getCollection,MbioDataset,character-method
setMethod("getCollection", "MbioDataset", function(object, collectionName = character(0), format = c("AbundanceData", "phyloseq", "Collection"), continuousMetadataOnly = c(FALSE, TRUE)) {
    format <- veupathUtils::matchArg(format)
    continuousMetadataOnly <- veupathUtils::matchArg(continuousMetadataOnly)
    
    if (length(collectionName) == 0) {
        stop("Must specify a collection name")
    }
    
    if (!collectionName %in% getCollectionNames(object)) {
        stop(sprintf("Collection '%s' does not exist", collectionName))
    }

    collection <- object@collections[[which(getCollectionNames(object) == collectionName)]]
    if (format == "Collection") {
        return(collection)
    }

    collectionDT <- data.table::setDT(collection@data)
    collectionIdColumns <- c(collection@recordIdColumn, collection@ancestorIdColumns)

    # remove IRI from collection column names. strip everything after and including the last square bracket, and remove trailing spaces
    rawNames <- names(collectionDT)
    names(collectionDT)[! names(collectionDT) %in% collectionIdColumns] <- sub("\\s*\\[([^\\[]*)$", "", names(collectionDT)[! names(collectionDT) %in% collectionIdColumns])

    # provide as much taxonomic resolution as possible for vague cases
    unclassifiedIndexes <- which(grepl('Incertae Sedis', names(collectionDT), fixed=TRUE))
    names(collectionDT)[unclassifiedIndexes] <- rawNames[unclassifiedIndexes]

    if (!!length(object@metadata@data)) {

        # need to be sure sample metadata contains only the relevant rows, actually having assay data
        # also need to make sure it has the assay record id column
        sampleMetadataDT <- data.table::setDT(merge(
            object@metadata@data, 
            collectionDT[, collectionIdColumns, with = FALSE], 
            by = c(object@metadata@ancestorIdColumns, object@metadata@recordIdColumn),
            all.y = TRUE
        ))

        # if we only want continuous metadata, only keep numeric columns
        # means we lose dates, but i think thats ok for now
        if (continuousMetadataOnly) {
            metadataColNames <- names(sampleMetadataDT)
            numericColumns <- metadataColNames[which(sapply(sampleMetadataDT,is.numeric))]
            sampleMetadataDT <- sampleMetadataDT[, unique(c(collectionIdColumns, numericColumns)), with=FALSE]
        }

        # also need to make sure they are in the same order
        data.table::setorderv(sampleMetadataDT, cols=collection@recordIdColumn)
        data.table::setorderv(collectionDT, cols=collection@recordIdColumn)

        sampleMetadata <- new("SampleMetadata",
            data = sampleMetadataDT,
            recordIdColumn = collection@recordIdColumn,
            ancestorIdColumns = collection@ancestorIdColumns
        )

    } else {
        sampleMetadataDT <- data.table::data.table()
        sampleMetadata <- SampleMetadata()
    }
    
    if (format == "AbundanceData" || format == "DESeqDataSet") {

        collectionDataDT <- collectionDT[, -collectionIdColumns, with = FALSE]
        if (all(collectionDataDT == round(collectionDataDT), na.rm = TRUE)) {
            veupathUtils::logWithTime("Integer values detected. Converting collection to AbsoluteAbundanceData", verbose = TRUE)
            abundanceData <- AbsoluteAbundanceData(
                name = collection@name,
                data = collectionDT, 
                sampleMetadata = sampleMetadata, 
                recordIdColumn = collection@recordIdColumn,
                ancestorIdColumns = collection@ancestorIdColumns
            )
        } else {
            abundanceData <- AbundanceData(
                name = collection@name,
                data = collectionDT, 
                sampleMetadata = sampleMetadata, 
                recordIdColumn = collection@recordIdColumn,
                ancestorIdColumns = collection@ancestorIdColumns
            )
        }
    } else if (format == "phyloseq") {
        .require_package("phyloseq")

        sampleNames <- collectionDT[[collection@recordIdColumn]]
        keepCols <- names(collectionDT)[! names(collectionDT) %in% collectionIdColumns]
        taxaNames <- names(collectionDT[, keepCols, with = FALSE])
        otu <- t(collectionDT[, keepCols, with = FALSE])
        names(otu) <- sampleNames
        rownames(otu) <- taxaNames

        tax <- data.frame(taxonomy = rownames(otu))
        rownames(tax) <- tax$taxonomy

        samples <- sampleMetadataDT

        if (nrow(samples) != 0) {
            rownames(samples) <- sampleNames
            samples <- samples[, !collectionIdColumns, with = FALSE]

            abundanceData <- phyloseq::phyloseq(
                phyloseq::otu_table(as.matrix(otu), taxa_are_rows = TRUE),
                phyloseq::sample_data(samples),
                phyloseq::tax_table(as.matrix(tax))
            )
        } else {
            abundanceData <- phyloseq::phyloseq(
                phyloseq::otu_table(as.matrix(otu), taxa_are_rows = TRUE),
                phyloseq::tax_table(as.matrix(tax))
            )
        }
    } 

    return(abundanceData)
})

collectionVarNamesGeneric <- getGeneric("getCollectionVariableNames", "veupathUtils")
#' Get Microbiome Dataset Collection Variable Names
#' 
#' Get the variable names in a collection in the Microbiome Dataset.
#' 
#' @examples
#' variableNames <- getCollectionVariableNames(
#'	microbiomeData::DiabImmune, 
#'	"16S (V4) Genus (Relative taxonomic abundance analysis)"
#' )
#' @param object A Microbiome Dataset
#' @param collectionName The name of the collection to return the variable names for
#' @return a character vector of the variable names in the requested collection
#' @export
#' @importFrom veupathUtils getCollectionVariableNames
setMethod(collectionVarNamesGeneric, "MbioDataset", function(object, collectionName) {
    return(veupathUtils::getCollectionVariableNames(getCollection(object, collectionName)))
})

#' Get Microbiome Dataset Variables
#' 
#' Get the variables in the Microbiome Dataset by their names.
#' The requested variables could belong to any collection or
#' to the metadata. The returned data.table will contain the 
#' requested variables as columns and any appropriate identifiers. 
#' If one of the requested variables cannot be returned, a warning
#' will be printed.
#' 
#' @examples
#' getCollectionVariableNames(
#'	microbiomeData::DiabImmune, 
#'	"16S (V4) Genus (Relative taxonomic abundance analysis)"
#' )
#'
#' getMetadataVariableNames(microbiomeData::DiabImmune)
#' variablesDT <- getVariables(
#'      microbiomeData::DiabImmune, 
#'      list("metadata" = c("age_months", "sex"),
#'           "16S (V4) Genus (Relative taxonomic abundance analysis)" = "Bacteroides", 
#'           "Shotgun metagenomics Metagenome enzyme pathway abundance data" = "ANAGLYCOLYSIS-PWY: glycolysis III (from glucose)"
#'      )
#' )
#' @param object A Microbiome Dataset
#' @param variables The names of the variables to return. This should be a named list
#' where the names are collection names and the values are variable names for that collection.
#' For the case of metadata variables, the name should be "metadata".
#' @return a data.table of the requested variables
#' @rdname getVariables
#' @export
setGeneric("getVariables", function(object, variables) standardGeneric("getVariables"), signature = "object")

#' @rdname getVariables
#' @aliases getVariables,MbioDataset,character-method
setMethod("getVariables", "MbioDataset", function(object, variables) {

    if (!is.list(variables)) {
        stop("variables argument must be a list")
    }
    if (is.null(names(variables))) {
        stop("variables argument must be a named list")
    }

    ## identify variables w identical names early
    flattenedVars <- unlist(variables)
    dups <- unname(flattenedVars[duplicated(flattenedVars)])
    collectionsWithDups <- lapply(variables, function(x) {dups %in% x})
    collectionsWithDupsIndexes <- unname(which(collectionsWithDups == TRUE))

    fetchCollectionVariables <- function(collectionIndex) {
        variableNames <- variables[[collectionIndex]]
        collectionName <- names(variables)[collectionIndex]

        if (collectionName == "metadata") {
            return(getSampleMetadata(object, metadataVariables = variableNames))
        }

        if (!collectionName %in% getCollectionNames(object)) {
            stop(sprintf("Collection '%s' does not exist", collectionName))
        }

        if (any(variableNames %in% getCollectionVariableNames(object, collectionName))) {
            collection <- getCollection(object, collectionName)
            presentVars <- variableNames[variableNames %in% getCollectionVariableNames(collection)]
            if (veupathUtils::isOneToManyWithAncestor(collection)) {
                warning("Unable to return the following variables: ", presentVars)
                return(data.table::data.table())
            }
            dt <- veupathUtils::getCollectionData(collection, presentVars)
            if (collectionIndex %in% collectionsWithDupsIndexes) {
                ## rename variables to prepend the collection name
                names(dt)[names(dt) %in% presentVars] <- paste(collectionName, presentVars)
            }
            return(dt)
        } else {
            stop(sprintf("Collection '%s' does not contain the following variables: %s", collectionName, paste(variableNames, collapse = ", ")))
        }
    }

    ## this kind of assumes that metadata are always on ancestor entities of assays
    ## this will break for user data, when we get there
    mergeCols <- getSampleMetadataIdColumns(object)

    if (length(variables) == 0) {
        return(data.table::data.table())
    }

    collectionVarDTs <- lapply(1:length(variables), fetchCollectionVariables)
    names(collectionVarDTs) <- names(variables)
    collectionVarDT <- purrr::reduce(collectionVarDTs, customMerge, mergeCols = mergeCols)

    return(collectionVarDT)
})

## a helper that merges two collections of variables
## if either input is empty, returns the other
## use this w some caution. It is barely a general
## purpose function, and isnt really tested.
customMerge <- function(x, y, mergeCols = NULL) {
    if (!inherits(x, "data.table")) {
        stop("Argument 'x' must be a data.table")
    } else if (!inherits(y, "data.table")) {
        stop("Argument 'y' must be a data.table")
    }

    if (!length(x)) {
        return(y)
    } else if (!length(y)) {
        return(x)
    } else {
        if (is.null(mergeCols)) {
            return(merge(x, y))
        }
        return(merge(x, y, by = mergeCols))
    }
}
