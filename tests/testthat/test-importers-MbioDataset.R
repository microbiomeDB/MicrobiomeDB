test_that("buildCollectionFromTreeSE works", {
    # should return a Collection, and validate the normalization is accurate

    # use simple mock data for this
    # treeSE should have matrix class for assay data and data.frame class for row data
    assay_data <- rbind(rep(0, 4), matrix(1:20, nrow = 5))
    colnames(assay_data) <- paste0("sample", 1:4)
    rownames(assay_data) <- paste("entity", seq_len(6), sep = "")

    row_data <- data.frame(Kingdom = "A",
                        Phylum = rep(c("B1", "B2"), c(2, 4)),
                        Class = rep(c("C1", "C2", "C3"), each = 2),
                        OTU = paste0("D", 1:6),
                        row.names = rownames(assay_data),
                        stringsAsFactors = FALSE)

    # no normalization
    collectionRaw <- buildCollectionFromTreeSE(
        collectionName = list(assayDataName = "test", rowDataColumnName = "OTU"),
        rowData = row_data,
        assayData = assay_data,
        normalizationMethod = "none",
        verbose = TRUE
    )

    expect_equal(inherits(collectionRaw, "Collection"), TRUE)
    expect_equal(collectionRaw@name, "test: OTU")
    expect_equal(length(collectionRaw@data), 7) # 6 OTUs + 1 recordIDs

    # TSS normalized
    collectionNormalized <- buildCollectionFromTreeSE(
        collectionName = list(assayDataName = "test", rowDataColumnName = "OTU"),
        rowData = row_data,
        assayData = assay_data,
        normalizationMethod = "TSS",
        verbose = TRUE
    )

    expect_equal(inherits(collectionNormalized, "Collection"), TRUE)
    expect_equal(collectionNormalized@name, "test: OTU (TSS normalized)")
    expect_equal(length(collectionNormalized@data), 7)
    expect_equal(all(rowSums(collectionNormalized@data[, -"recordIDs"]) == 1), TRUE)

    # try with Class, make sure its aggregating OTU to the Class level
    collectionClass <- buildCollectionFromTreeSE(
        collectionName = list(assayDataName = "test", rowDataColumnName = "Class"),
        rowData = row_data,
        assayData = assay_data,
        normalizationMethod = "TSS",
        verbose = TRUE
    )

    expect_equal(inherits(collectionClass, "Collection"), TRUE)
    expect_equal(collectionClass@name, "test: Class (TSS normalized)")
    expect_equal(length(collectionClass@data), 4)
    expect_equal(all(rowSums(collectionClass@data[, -"recordIDs"]) == 1), TRUE)
})

test_that("we can get an MbioDataset from a TreeSummarizedExperiment", {
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns

    # no normalization, with raw values
    mbioDataset <- importTreeSummarizedExperiment(tse, normalizationMethod = "none", keepRawValues = TRUE, verbose = TRUE)
    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE)
    expect_equal("counts: Genus" %in% getCollectionNames(mbioDataset), TRUE)
    expect_equal("counts: Genus (TSS normalized)" %in% getCollectionNames(mbioDataset), FALSE)
    expect_equal(length(getCollection(mbioDataset, "counts: Genus")@data) > 1, TRUE)

    # TSS normalization, no raw values
    mbioDataset <- importTreeSummarizedExperiment(tse, normalizationMethod = "TSS", keepRawValues = FALSE, verbose = TRUE)
    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE)
    expect_equal("counts: Genus (TSS normalized)" %in% getCollectionNames(mbioDataset), TRUE)
    expect_equal("counts: Genus" %in% getCollectionNames(mbioDataset), FALSE)
    expect_equal(length(getCollection(mbioDataset, "counts: Genus (TSS normalized)")@data) > 1, TRUE)

    # TSS normalization, with raw values
    mbioDataset <- importTreeSummarizedExperiment(tse, normalizationMethod = "TSS", keepRawValues = TRUE, verbose = TRUE)
    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE)
    expect_equal("counts: Genus (TSS normalized)" %in% getCollectionNames(mbioDataset), TRUE)
    expect_equal("counts: Genus" %in% getCollectionNames(mbioDataset), TRUE)
    expect_equal(length(getCollection(mbioDataset, "counts: Genus (TSS normalized)")@data) > 1, TRUE)
    expect_equal(length(getCollection(mbioDataset, "counts: Genus")@data) > 1, TRUE)

    # no normalization, no raw values (should err)
    expect_error(importTreeSummarizedExperiment(tse, normalizationMethod = "none", keepRawValues = FALSE, verbose = TRUE))

    # not a summarized experiment
    expect_error(importTreeSummarizedExperiment(data.frame(), normalizationMethod = "none", keepRawValues = FALSE, verbose = TRUE))

})

# most of the important testing is in mia, so 
# we can test we get the right class back and its populated..
# TODO add instructions to readme for installing mia, bc of the system level dep
test_that("the humann miaverse wrapper works", {
    skip_if_not_installed("mia")

    file_path <- system.file("extdata", "humann_output.tsv", package = "mia")

    mbioDataset <- importHUMAnN(normalizationMethod = "none", keepRawValues = TRUE, verbose = TRUE, file_path)

    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE)
    expect_equal(length(getCollectionNames(mbioDataset)) > 0, TRUE)
    aCollectionName <- getCollectionNames(mbioDataset)[1]
    aCollection <- getCollection(mbioDataset, aCollectionName)
    expect_equal(inherits(aCollection, "Collection"), TRUE)
    expect_equal(length(aCollection@data) > 0, TRUE)

})

test_that("the mothur miaverse wrapper works", {
    skip_if_not_installed("mia")
    
    counts <- system.file("extdata", "mothur_example.shared", package = "mia")
    taxa <- system.file("extdata", "mothur_example.cons.taxonomy", package = "mia")
    taxa2 <- system.file("extdata", "mothur_example.taxonomy", package = "mia")
    meta <- system.file("extdata", "mothur_example.design", package = "mia")

    mbioDataset <- importMothur(normalizationMethod = "none", keepRawValues = TRUE, verbose = TRUE, counts)

    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE)
    expect_equal(length(getCollectionNames(mbioDataset)) > 0, FALSE) # no real data, needs taxa file too

    mbioDataset <- importMothur(normalizationMethod = "none", keepRawValues = TRUE, verbose = TRUE, counts, taxa)

    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE) # this one should work but not be interesting without meta
    expect_equal(length(getCollectionNames(mbioDataset)) > 0, TRUE)
    aCollectionName <- getCollectionNames(mbioDataset)[1]
    aCollection <- getCollection(mbioDataset, aCollectionName)
    expect_equal(inherits(aCollection, "Collection"), TRUE)
    expect_equal(length(aCollection@data) > 0, TRUE)

    mbioDataset <- importMothur(normalizationMethod = "none", keepRawValues = TRUE, verbose = TRUE, counts, taxa2)

    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE) # this one should work but not be interesting without meta
    expect_equal(length(getCollectionNames(mbioDataset)) > 0, TRUE)
    aCollectionName <- getCollectionNames(mbioDataset)[1]
    aCollection <- getCollection(mbioDataset, aCollectionName)
    expect_equal(inherits(aCollection, "Collection"), TRUE)
    expect_equal(length(aCollection@data) > 0, TRUE)

    mbioDataset <- importMothur(normalizationMethod = "none", keepRawValues = TRUE, verbose = TRUE, counts, taxa, meta)

    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE) # this is the one we want really..
    expect_equal(length(getCollectionNames(mbioDataset)) > 0, TRUE)
    aCollectionName <- getCollectionNames(mbioDataset)[1]
    aCollection <- getCollection(mbioDataset, aCollectionName)
    expect_equal(inherits(aCollection, "Collection"), TRUE)
    expect_equal(length(aCollection@data) > 0, TRUE)

})

test_that("the qiime2 miaverse wrapper works", {
    skip_if_not_installed("mia")
    
    featureTableFile <- system.file("extdata", "table.qza", package = "mia")
    taxonomyTableFile <- system.file("extdata", "taxonomy.qza", package = "mia")

    mbioDataset <- importQIIME2(normalizationMethod = "none", keepRawValues = TRUE, verbose = TRUE, featureTableFile, taxonomyTableFile)

    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE)
    expect_equal(length(getCollectionNames(mbioDataset)) > 0, TRUE)
    aCollectionName <- getCollectionNames(mbioDataset)[1]
    aCollection <- getCollection(mbioDataset, aCollectionName)
    expect_equal(inherits(aCollection, "Collection"), TRUE)
    expect_equal(length(aCollection@data) > 0, TRUE)

})

test_that("the biom miaverse wrapper works", {
    skip_if_not_installed("mia")
    skip_if_not_installed("biomformat")
    
    rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom",
                                    package = "biomformat")

    mbioDataset <- importBIOM(normalizationMethod = "none", keepRawValues = TRUE, verbose = TRUE, rich_dense_file)

    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE)
    expect_equal(length(getCollectionNames(mbioDataset)) > 0, TRUE)
    aCollectionName <- getCollectionNames(mbioDataset)[1]
    aCollection <- getCollection(mbioDataset, aCollectionName)
    expect_equal(inherits(aCollection, "Collection"), TRUE)
    expect_equal(length(aCollection@data) > 0, TRUE)

    rich_dense_biom = biomformat::read_biom(rich_dense_file)

    mbioDataset <- importBIOM(normalizationMethod = "none", keepRawValues = TRUE, verbose = TRUE, rich_dense_biom)

    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE)
    expect_equal(length(getCollectionNames(mbioDataset)) > 0, TRUE)
    aCollectionName <- getCollectionNames(mbioDataset)[1]
    aCollection <- getCollection(mbioDataset, aCollectionName)
    expect_equal(inherits(aCollection, "Collection"), TRUE)
    expect_equal(length(aCollection@data) > 0, TRUE)
})

test_that("the dada2 miaverse wrapper works", {
    skip_if_not_installed("mia")
    skip_if_not_installed("dada2")

    fnF <- system.file("extdata", "sam1F.fastq.gz", package="dada2")
    fnR = system.file("extdata", "sam1R.fastq.gz", package="dada2")
    dadaF <- dada2::dada(fnF, selfConsist=TRUE)
    dadaR <- dada2::dada(fnR, selfConsist=TRUE)

    mbioDataset <- importDADA2(normalizationMethod = "none", keepRawValues = TRUE, verbose = TRUE, dadaF, fnF, dadaR, fnR)

    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE)
    expect_equal(length(getCollectionNames(mbioDataset)) > 0, FALSE) # dummy daddy

})

test_that("the phyloseq miaverse wrapper works", {
    skip_if_not_installed("mia")
    skip_if_not_installed("phyloseq")

    data(GlobalPatterns, package="phyloseq")

    mbioDataset <- importPhyloseq(normalizationMethod = "none", keepRawValues = TRUE, verbose = TRUE, GlobalPatterns)

    expect_equal(inherits(mbioDataset, "MbioDataset"), TRUE)
    expect_equal(length(getCollectionNames(mbioDataset)) > 0, TRUE)
    aCollectionName <- getCollectionNames(mbioDataset)[1]
    aCollection <- getCollection(mbioDataset, aCollectionName)
    expect_equal(inherits(aCollection, "Collection"), TRUE)
    expect_equal(length(aCollection@data) > 0, TRUE)
})