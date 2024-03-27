collectionNamesGeneric <- getGeneric("getCollectionNames", "veupathUtils")
#' Get Names of Collections
#' 
#' Get the names of the collections in a MbioDataset object
#' @param object An MbioDataset
#' @return A character vector of collection names
#' @export
setMethod(collectionNamesGeneric, "MbioDataset", function(object) return(unname(getCollectionNames(object@collections))))

metadataVarNamesGeneric <- getGeneric("getMetadataVariableNames", "veupathUtils")
#' Get Variable Names of Metadata
#' 
#' Get the names of the metadata variables in an MbioDataset.
#' @param object An MbioDataset
#' @return a character vector of metadata variable names
#' @export
setMethod(metadataVarNamesGeneric, "MbioDataset", function(object) return(names(object@metadata@data)))

sampleMetadataGeneric <- getGeneric("getSampleMetadata", "veupathUtils")
#' Get data.table of sample metadata from MbioDataset
#'
#' Returns a data.table of sample metadata
#' 
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
#' @param object A Microbiome Dataset
#' @param collectionName The name of the collection to return
#' @param format The format of the collection to return. Currently supported options are "AbundanceData", "phyloseq", and "Collection".
#' @param continuousMetadataOnly If TRUE, only continuous metadata will be returned. If FALSE, all metadata will be returned.
#' @return An AbundanceData, phyloseq, or Collection object representing the collection and any associated study metadata
#' @importFrom phyloseq phyloseq
#' @include class-AbundanceData.R
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
    names(collectionDT)[names(collectionDT) == 'Incertae Sedis'] <- rawNames[names(collectionDT) == 'Incertae Sedis']

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
    
    if (format == "AbundanceData") {

        abundanceData <- AbundanceData(
            name = collection@name,
            data = collectionDT, 
            sampleMetadata = sampleMetadata, 
            recordIdColumn = collection@recordIdColumn,
            ancestorIdColumns = collection@ancestorIdColumns
        )

    } else if (format == "phyloseq") {

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