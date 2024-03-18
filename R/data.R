## We want to be rather specific about what data and methods from the data and computes packages
## we make available to users when they attach this package. To that end, Im explicitly importing
## and re-exporting some things, rather than addings these two packages in Depends in DESCRIPTION


#' List Curated Datasets
#'
#' This function lists curated datasets from MicrobiomeDB.org which are available in this package.
#' @export
getCuratedDatasetNames <- function() {

	## TODO try this instead, though idk if attaching the monster is a good idea
	# also better to put it in the data package, import and reexport here
	#.filename <- system.file("R", "sysdata.rda", package = "microbiomeData")
	#attached_filename <- paste0("file:", .filename, "")
    #suppressMessages(do.call("attach", list(what = .filename, name = attached_filename)))
    #on.exit(eval(substitute(detach(name), list(name = attached_filename))))
    #return(ls(envir = as.environment(attached_filename)))
	
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

#' @importFrom microbiomeData GEMS1
#' @export
microbiomeData::GEMS1

#' @importFrom microbiomeData DailyBaby
#' @export
microbiomeData::DailyBaby

#' @importFrom microbiomeData HMP_V1V3
#' @export
microbiomeData::HMP_V1V3

#' @importFrom microbiomeData HMP_V3V5
#' @export
microbiomeData::HMP_V3V5

#' @importFrom microbiomeData Anopheles_albimanus
#' @export
microbiomeData::Anopheles_albimanus

#' @importFrom microbiomeData ECAM
#' @export
microbiomeData::ECAM

#' @importFrom microbiomeData MORDOR
#' @export
microbiomeData::MORDOR

#' @importFrom microbiomeData EcoCF
#' @export
microbiomeData::EcoCF

#' @importFrom microbiomeData Leishmaniasis
#' @export
microbiomeData::Leishmaniasis

#' @importFrom microbiomeData MALED_2yr
#' @export
microbiomeData::MALED_2yr

#' @importFrom microbiomeData MALED_diarrhea
#' @export
microbiomeData::MALED_diarrhea

#' @importFrom microbiomeData Malaysia_helminth
#' @export
microbiomeData::Malaysia_helminth

#' @importFrom microbiomeData PIH_Uganda
#' @export
microbiomeData::PIH_Uganda

#' @importFrom microbiomeData PretermInfantResistome1
#' @export
microbiomeData::PretermInfantResistome1