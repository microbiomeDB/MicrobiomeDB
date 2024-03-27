check_mbio_dataset <- function(object) {
    errors <- character()

    # check that at least some ancestorIdColumns are shared between collections and sampleMetadata
    sampleMetadataAncestorIds <- object@metadata@ancestorIdColumns
    if (!!length(sampleMetadataAncestorIds)) {
        if (!all(sapply(object@collections, function(x) any(x@ancestorIdColumns %in% sampleMetadataAncestorIds)))) {
            msg <- "at least one ancestorIdColumn must be shared between collections and sampleMetadata"
            errors <- c(errors, msg)
        }
    }
    
    if (length(errors) == 0) {
        return(TRUE)
    } else {
        return(errors)
    }
}

#' MicrobiomeDB Dataset
#' 
#' This class represents a MicrobiomeDB dataset.
#' @name MbioDataset-class
#' @rdname MbioDataset-class
#' @importFrom veupathUtils SampleMetadata
setClass("MbioDataset", 
    slots = c(
        collections = "Collections",
        metadata = "SampleMetadata"
    ),
    validity = check_mbio_dataset
)