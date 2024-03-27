

sampleMetadataBuilder <- function(dataSource) {
    dt <- getDataFromSource(dataSource, keepIdsAndNumbersOnly=FALSE, cleanColumnNames=TRUE)
    dataColNames <- names(dt)
    recordIdColumn <- findRecordIdColumn(dataColNames)
    ancestorIdColumns <- findAncestorIdColumns(dataColNames)

    sampleMetadata <- new("SampleMetadata",
        data=dt,
        recordIdColumn = recordIdColumn,
        ancestorIdColumns = ancestorIdColumns
    )

    return(sampleMetadata)
}

mergeSampleMetadata <- function(x, y) {
    uniqueAncestorIdColumns <- unique(c(x@ancestorIdColumns, y@ancestorIdColumns))
    recordIdColumn <- ifelse(x@recordIdColumn %in% uniqueAncestorIdColumns, y@recordIdColumn, x@recordIdColumn)
    data <- merge(x@data, y@data, by = uniqueAncestorIdColumns, all = TRUE)

    sampleMetadata <- new("SampleMetadata",
        data = data,
        recordIdColumn = recordIdColumn,
        ancestorIdColumns = uniqueAncestorIdColumns
    )

    return(sampleMetadata)
}

#' @importFrom purrr reduce
sampleMetadataFromDataSources <- function(dataSources) {
    sampleMetataList <- lapply(dataSources, sampleMetadataBuilder)
    sampleMetadata <- purrr::reduce(sampleMetataList, mergeSampleMetadata)

    return(sampleMetadata)
}


##### I have never hated S4 so much... even ai assisted this was excruciating #####

#' Create a Microbiome Dataset
#' 
#' This is a constructor for the MbioDataset class. It creates a MbioDataset containing
#' a list of Collections and a SampleMetadata object.
#' @param collections A list of Collection objects, a data.frame containing multiple collections,
#'  or a character vector containing one or more file path(s)
#' @param metadata A SampleMetadata object, a data.frame containing sample metadata,
#' or a list of file path(s)
#' @param ontology An data.frame containing the ontology of the dataset, or a character vector
#'  containing a file path to a data.frame
#' @export
#' @rdname MbioDataset-methods
setGeneric("MbioDataset", function(collections, metadata, ontology) standardGeneric("MbioDataset"))

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,missing,missing,missing-method
setMethod("MbioDataset", signature("missing", "missing", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset")
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,SampleMetadata,missing-method
setMethod("MbioDataset", signature("Collections", "SampleMetadata", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = collections, metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,SampleMetadata,data.frame-method
setMethod("MbioDataset", signature("Collections", "SampleMetadata", "data.frame"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collections object is provided.")
    new("MbioDataset", collections = collections, metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,SampleMetadata,character-method
setMethod("MbioDataset", signature("Collections", "SampleMetadata", "character"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collections object is provided.")
    new("MbioDataset", collections = collections, metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,data.frame,missing-method
setMethod("MbioDataset", signature("Collections", "data.frame", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = collections, metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,data.frame,data.frame-method
setMethod("MbioDataset", signature("Collections", "data.frame", "data.frame"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collections object is provided.")
    new("MbioDataset", collections = collections, metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,data.frame,character-method
setMethod("MbioDataset", signature("Collections", "data.frame", "character"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collections object is provided.")
    new("MbioDataset", collections = collections, metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,list,missing-method
setMethod("MbioDataset", signature("Collections", "list", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = collections, metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,list,data.frame-method
setMethod("MbioDataset", signature("Collections", "list", "data.frame"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collections object is provided.")
    new("MbioDataset", collections = collections, metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,list,character-method
setMethod("MbioDataset", signature("Collections", "list", "character"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collections object is provided.")
    new("MbioDataset", collections = collections, metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,missing,missing-method
setMethod("MbioDataset", signature("Collections", "missing", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = collections)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,missing,data.frame-method
setMethod("MbioDataset", signature("Collections", "missing", "data.frame"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collections object is provided.")
    new("MbioDataset", collections = collections, metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,missing,character-method
setMethod("MbioDataset", signature("Collections", "missing", "character"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collections object is provided.")
    new("MbioDataset", collections = collections, metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,character,missing-method
setMethod("MbioDataset", signature("Collections", "character", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = collections, metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,character,data.frame-method
setMethod("MbioDataset", signature("Collections", "character", "data.frame"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collections object is provided.")
    new("MbioDataset", collections = collections, metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collections,character,character-method
setMethod("MbioDataset", signature("Collections", "character", "character"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collections object is provided.")
    new("MbioDataset", collections = collections, metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,SampleMetadata,missing-method
setMethod("MbioDataset", signature("Collection", "SampleMetadata", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,SampleMetadata,data.frame-method
setMethod("MbioDataset", signature("Collection", "SampleMetadata", "data.frame"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collection object is provided.")
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,SampleMetadata,character-method
setMethod("MbioDataset", signature("Collection", "SampleMetadata", "character"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collection object is provided.")
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,data.frame,missing-method
setMethod("MbioDataset", signature("Collection", "data.frame", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,data.frame,data.frame-method
setMethod("MbioDataset", signature("Collection", "data.frame", "data.frame"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collection object is provided.")
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,data.frame,character-method
setMethod("MbioDataset", signature("Collection", "data.frame", "character"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collection object is provided.")
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,list,missing-method
setMethod("MbioDataset", signature("Collection", "list", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,list,data.frame-method
setMethod("MbioDataset", signature("Collection", "list", "data.frame"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collection object is provided.")
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,list,character-method
setMethod("MbioDataset", signature("Collection", "list", "character"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collection object is provided.")
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,missing,missing-method
setMethod("MbioDataset", signature("Collection", "missing", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,missing,data.frame-method
setMethod("MbioDataset", signature("Collection", "missing", "data.frame"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collection object is provided.")
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,missing,character-method
setMethod("MbioDataset", signature("Collection", "missing", "character"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collection object is provided.")
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,character,missing-method
setMethod("MbioDataset", signature("Collection", "character", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,character,data.frame-method
setMethod("MbioDataset", signature("Collection", "character", "data.frame"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collection object is provided.")
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,Collection,character,character-method
setMethod("MbioDataset", signature("Collection", "character", "character"), function(collections, metadata, ontology) {
    warning("Ontology specified but not used when a Collection object is provided.")
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,SampleMetadata,missing-method
setMethod("MbioDataset", signature("list", "SampleMetadata", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,SampleMetadata,data.frame-method
setMethod("MbioDataset", signature("list", "SampleMetadata", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,SampleMetadata,character-method
setMethod("MbioDataset", signature("list", "SampleMetadata", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,data.frame,missing-method
setMethod("MbioDataset", signature("list", "data.frame", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,data.frame,data.frame-method
setMethod("MbioDataset", signature("list", "data.frame", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,data.frame,character-method
setMethod("MbioDataset", signature("list", "data.frame", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,missing,missing-method
setMethod("MbioDataset", signature("list", "missing", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,missing,data.frame-method
setMethod("MbioDataset", signature("list", "missing", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,missing,character-method
setMethod("MbioDataset", signature("list", "missing", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,list,missing-method
setMethod("MbioDataset", signature("list", "list", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,list,data.frame-method
setMethod("MbioDataset", signature("list", "list", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,list,character-method
setMethod("MbioDataset", signature("list", "list", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,character,missing-method
setMethod("MbioDataset", signature("list", "character", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,character,data.frame-method
setMethod("MbioDataset", signature("list", "character", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,list,character,character-method
setMethod("MbioDataset", signature("list", "character", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,SampleMetadata,missing-method
setMethod("MbioDataset", signature("data.frame", "SampleMetadata", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,SampleMetadata,data.frame-method
setMethod("MbioDataset", signature("data.frame", "SampleMetadata", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,SampleMetadata,character-method
setMethod("MbioDataset", signature("data.frame", "SampleMetadata", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,list,missing-method
setMethod("MbioDataset", signature("data.frame", "list", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,list,data.frame-method
setMethod("MbioDataset", signature("data.frame", "list", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,list,character-method
setMethod("MbioDataset", signature("data.frame", "list", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,missing,missing-method
setMethod("MbioDataset", signature("data.frame", "missing", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,missing,data.frame-method
setMethod("MbioDataset", signature("data.frame", "missing", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,missing,character-method
setMethod("MbioDataset", signature("data.frame", "missing", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,data.frame,missing-method
setMethod("MbioDataset", signature("data.frame", "data.frame", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,data.frame,data.frame-method
setMethod("MbioDataset", signature("data.frame", "data.frame", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,data.frame,character-method
setMethod("MbioDataset", signature("data.frame", "data.frame", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,character,missing-method
setMethod("MbioDataset", signature("data.frame", "character", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,character,data.frame-method
setMethod("MbioDataset", signature("data.frame", "character", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,data.frame,character,character-method
setMethod("MbioDataset", signature("data.frame", "character", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,character,missing-method
setMethod("MbioDataset", signature("character", "character", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,character,data.frame-method
setMethod("MbioDataset", signature("character", "character", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,character,character-method
setMethod("MbioDataset", signature("character", "character", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = sampleMetadataFromDataSources(list(metadata)))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,list,missing-method
setMethod("MbioDataset", signature("character", "list", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,list,data.frame-method
setMethod("MbioDataset", signature("character", "list", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,list,character-method
setMethod("MbioDataset", signature("character", "list", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = sampleMetadataFromDataSources(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,missing,missing-method
setMethod("MbioDataset", signature("character", "missing", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,missing,data.frame-method
setMethod("MbioDataset", signature("character", "missing", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,missing,character-method
setMethod("MbioDataset", signature("character", "missing", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,data.frame,missing-method
setMethod("MbioDataset", signature("character", "data.frame", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,data.frame,data.frame-method
setMethod("MbioDataset", signature("character", "data.frame", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,data.frame,character-method
setMethod("MbioDataset", signature("character", "data.frame", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = sampleMetadataBuilder(metadata))
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,SampleMetadata,missing-method
setMethod("MbioDataset", signature("character", "SampleMetadata", "missing"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections), metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,SampleMetadata,data.frame-method
setMethod("MbioDataset", signature("character", "SampleMetadata", "data.frame"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, ontology), metadata = metadata)
})

#' @rdname MbioDataset-methods
#' @aliases MbioDataset,character,SampleMetadata,character-method
setMethod("MbioDataset", signature("character", "SampleMetadata", "character"), function(collections, metadata, ontology) {
    new("MbioDataset", collections = veupathUtils::Collections(collections, data.table::fread(ontology)), metadata = metadata)
})