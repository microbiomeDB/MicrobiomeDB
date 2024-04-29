test_that("we can create a new MbioDataset", {
    # an empty one
    mbioDataset <- MbioDataset()
    expect_s4_class(mbioDataset, "MbioDataset")
    
    # a manually populated one
    mbioDataset <- MbioDataset(
        Collections(list(Collection("my collection", data.frame(entity.id = 1, entity.collection_x = 1, entity.collection_y = 2, ancestor.y = 1), "entity.id", "ancestor.y"))), 
        SampleMetadata()
    )
    expect_s4_class(mbioDataset, "MbioDataset")   

    # from a Collection object and SampleMetadata
    mbioDataset <- MbioDataset(
        Collection("my collection", data.frame(entity.id = 1, entity.collection_x = 1, entity.collection_y = 2, ancestor.y = 1), "entity.id", "ancestor.y"), 
        SampleMetadata()
    )
    expect_s4_class(mbioDataset, "MbioDataset")
    
    # from a file of collections and sample metadata
    dataFile1 <- test_path('testdata','DiabImmune/DiabImmune_entity_16SRRNAV4Assay.txt')
    metadataFile1 <- test_path('testdata','DiabImmune/DiabImmune_ParticipantRepeatedMeasure.txt')
    mbioDataset <- MbioDataset(dataFile1, metadataFile1)
    expect_s4_class(mbioDataset, "MbioDataset")

    # from a data frame of collections and sample metadata
    df <- data.table::fread(dataFile1)
    metadata <- data.table::fread(metadataFile1)
    mbioDataset <- MbioDataset(df, metadata)
    expect_s4_class(mbioDataset, "MbioDataset")

    # from a list of files of collections and sample metadata
    dataFile2 <- test_path('testdata','DiabImmune/DiabImmune_MetagenomicSequencingAssay.txt')
    metadataFile2 <- test_path('testdata','DiabImmune/DiabImmune_Participant.txt')
    metadataFile3 <- test_path('testdata','DiabImmune/DiabImmune_Sample.txt')
    mbioDataset <- MbioDataset(list(dataFile1, dataFile2), list(metadataFile1, metadataFile2, metadataFile3))
    expect_s4_class(mbioDataset, "MbioDataset")

    # from a list of collections and sample metadata data frames
    df2 <- data.table::fread(dataFile2)
    metadata2 <- data.table::fread(metadataFile2)
    metadata3 <- data.table::fread(metadataFile3)
    mbioDataset <- MbioDataset(list(df, df2), list(metadata, metadata2, metadata3))
    expect_s4_class(mbioDataset, "MbioDataset")
})

test_that("we can update collection names and get collections", {
    mbioDataset <- MbioDataset(
        list(
            Collection("my collection", 
                data.frame(entity.id = 1, entity.collection_x = 1, entity.collection_y = 2, ancestor.y = 1), "entity.id", "ancestor.y"), 
            Collection("my collection 2", 
                data.frame(entity.id = 1, entity.collection_x = .1, entity.collection_y = .2, ancestor.y = 1), "entity.id", "ancestor.y")
        ),
        SampleMetadata()
    )

    testDataset <- updateCollectionName(mbioDataset, "my collection", "My Collection")

    expect_equal(testDataset@collections[[1]]@name, "My Collection")
    expect_equal(getCollectionNames(testDataset)[[1]], "My Collection")

    testCollection <- getCollection(testDataset, "My Collection")
    expect_s4_class(testCollection, "AbsoluteAbundanceData")
    expect_equal(testCollection@data, data.table::data.table(entity.id = 1, entity.collection_x = 1, entity.collection_y = 2, ancestor.y = 1))
    expect_equal(testCollection@recordIdColumn, "entity.id")
    expect_equal(testCollection@ancestorIdColumns, "ancestor.y")

    testCollection <- getCollection(testDataset, "my collection 2")
    expect_s4_class(testCollection, "AbundanceData")
    expect_equal(testCollection@data, data.table::data.table(entity.id = 1, entity.collection_x = .1, entity.collection_y = .2, ancestor.y = 1))
    expect_equal(testCollection@recordIdColumn, "entity.id")
    expect_equal(testCollection@ancestorIdColumns, "ancestor.y")

    testCollection <- getCollection(testDataset, "My Collection", "Collection")
    expect_s4_class(testCollection, "Collection")

    testCollection <- getCollection(testDataset, "My Collection", "phyloseq")
    expect_s4_class(testCollection, "phyloseq")
})

test_that("we can get arbitrary variables", {
    dataFile1 <- test_path('testdata','DiabImmune/DiabImmune_entity_16SRRNAV4Assay.txt')
    metadataFile1 <- test_path('testdata','DiabImmune/DiabImmune_ParticipantRepeatedMeasure.txt')
    dataFile2 <- test_path('testdata','DiabImmune/DiabImmune_MetagenomicSequencingAssay.txt')
    metadataFile2 <- test_path('testdata','DiabImmune/DiabImmune_Participant.txt')
    metadataFile3 <- test_path('testdata','DiabImmune/DiabImmune_Sample.txt')
    ontologyFile <- test_path('testdata','DiabImmune/DiabImmune_OntologyMetadata.txt')
    mbioDataset <- MbioDataset(list(dataFile1, dataFile2), list(metadataFile2, metadataFile1, metadataFile3), ontologyFile)

    # try a sensible thing w vars on different 1:1 entities
    variablesDT <- getVariables(
        mbioDataset, 
        list("metadata" = c("age_months", "sex"),
            "16S (V4) Genus" = "Bacteroides",
            "Shotgun metagenomics Metagenome enzyme pathway abundance data" = "ANAGLYCOLYSIS-PWY: glycolysis III (from glucose)"
            )
    )
    # expect a data.table w four columns
    expect_s3_class(variablesDT, "data.table")
    expect_equal(length(variablesDT), 9) # 4 vars + 5 ids
    expect_equal(all(c("age_months", "sex", "Bacteroides", "ANAGLYCOLYSIS-PWY: glycolysis III (from glucose)") %in% names(variablesDT)), TRUE)
    expect_equal(nrow(variablesDT) > 0, TRUE)

    # try a var that doesnt exist
    expect_error(
        variablesDT <- getVariables(
            mbioDataset, 
            list("metadata" = c("age_months", "sex"),
                "16S (V4) Genus" = "Bacteroides",
                "Shotgun metagenomics Metagenome enzyme pathway abundance data" = "ANAGLYCOLYSIS-PWY: glycolysis III (from glucose)",
                "Shotgun metagenomics Genus" = "doesntexist"
                )
        )
    )
    

    # try a collection that doesnt exist
    expect_error(
        variablesDT <- getVariables(
            mbioDataset, 
            list("metadata" = c("age_months", "sex"),
                "16S (V4) Genus" = "Bacteroides",
                "doesntexist" = "ANAGLYCOLYSIS-PWY: glycolysis III (from glucose)"
                )
        )
    )
    
    # try the same named variable on two different collections
    variablesDT <- getVariables(
        mbioDataset, 
        list("metadata" = c("age_months", "sex"),
            "16S (V4) Genus" = "Bacteroides",
            "Shotgun metagenomics Genus" = "Bacteroides"
            )
    )

    expect_s3_class(variablesDT, "data.table")
    expect_equal(length(variablesDT), 9) # 4 vars + 5 ids
    expect_equal(all(c("age_months", "sex", "16S (V4) Genus Bacteroides", "Shotgun metagenomics Genus Bacteroides") %in% names(variablesDT)), TRUE)
    expect_equal(nrow(variablesDT) > 0, TRUE)

    # pass something other than a named list
    expect_error(variablesDT <- getVariables(mbioDataset, "16S (V4) Genus"))
    expect_error(variablesDT <- getVariables(mbioDataset, list("16S (V4) Genus)")))

    # find an ex where assays arent 1:1 w ancestors
})