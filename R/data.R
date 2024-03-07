## We want to be rather specific about what data and methods from the data and computes packages
## we make available to users when they attach this package. To that end, Im explicitly importing
## and re-exporting some things, rather than addings these two packages in Depends in DESCRIPTION

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


#' @export
DiabImmune <- function() {
	microbiomeData::DiabImmune
}

#' @export
FARMM <- function() {
	microbiomeData::FARMM
}

#' @export
Bangladesh <- function() {
	microbiomeData::Bangladesh
}

#' @export
HMP_WGS <- function() {
	microbiomeData::HMP_WGS
}

#' @export
BONUS <- function() {
	microbiomeData::BONUS
}

#' @export
NICU_NEC <- function() {
	microbiomeData::NICU_NEC
}