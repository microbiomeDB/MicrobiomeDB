## We want to be rather specific about what data and methods from the data and computes packages
## we make available to users when they attach this package. To that end, Im explicitly importing
## and re-exporting some things, rather than addings these two packages in Depends in DESCRIPTION

# TODO figure a way to not have to specify these manually.. 
# maybe write this in the data package and reexport? then itd have access to sysdata.rda directly
#' List Curated Datasets
#'
#' This function lists curated datasets from MicrobiomeDB.org which are available in this package.
#' @export
getCuratedDatasetNames <- function() {
	c(
	    'DiabImmune',
	    'FARMM',
	    'Bangladesh',
	    'HMP_WGS',
	    'BONUS',
	    'NICU_NEC'
	)
}

#' @importFrom microbiomeData DiabImmune
#' @export
microbiomeData::DiabImmune

#' @importFrom microbiomeData FARMM
#' @export
microbiomeData::FARMM

#' @importFrom microbiomeData Bangladesh
#' @export
microbiomeData::Bangladesh

#' @importFrom microbiomeData HMP_WGS
#' @export
microbiomeData::HMP_WGS

#' @importFrom microbiomeData BONUS
#' @export
microbiomeData::BONUS

#' @importFrom microbiomeData NICU_NEC
#' @export
microbiomeData::NICU_NEC