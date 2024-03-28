test_that("we can get compute results in different formats", {
    dataFile1 <- test_path('testdata','DiabImmune/DiabImmune_entity_16SRRNAV4Assay.txt')
    metadataFile1 <- test_path('testdata','DiabImmune/DiabImmune_ParticipantRepeatedMeasure.txt')
    dataFile2 <- test_path('testdata','DiabImmune/DiabImmune_MetagenomicSequencingAssay.txt')
    metadataFile2 <- test_path('testdata','DiabImmune/DiabImmune_Participant.txt')
    metadataFile3 <- test_path('testdata','DiabImmune/DiabImmune_Sample.txt')
    ontologyFile <- test_path('testdata','DiabImmune/DiabImmune_OntologyMetadata.txt')
    mbioDataset <- MbioDataset(list(dataFile1, dataFile2), list(metadataFile2, metadataFile1, metadataFile3), ontologyFile)

    genus <- getCollection(mbioDataset, "16S Genus", continuousMetadataOnly = TRUE)

    # make sure metadata dont contain IRIs
    expect_equal(all(grepl('[',names(genus@sampleMetadata@data),fixed=T)), FALSE)
    # and that that means diff abund works now
    comparatorVariable <- microbiomeComputations::Comparator(
                        variable = veupathUtils::VariableMetadata(
                            variableSpec = veupathUtils::VariableSpec(
                                variableId = 'delivery_mode',
                                entityId = ''
                            ),
                            dataShape = veupathUtils::DataShape(value="BINARY")
                        ),
                        groupA = veupathUtils::BinList(
                            S4Vectors::SimpleList(
                                c(veupathUtils::Bin(
                                    binLabel="Vaginal"
                                ))
                            )
                        ),
                        groupB = veupathUtils::BinList(
                            S4Vectors::SimpleList(
                                c(veupathUtils::Bin(
                                    binLabel="Cesarean"
                                ))
                            )
                        )
    )
    diffAbundOutput <- microbiomeComputations::differentialAbundance(getCollection(mbioDataset, "16S Genus"), comparatorVariable, method='Maaslin2', verbose=FALSE)
    expect_equal(inherits(diffAbundOutput, "ComputeResult"), TRUE)

    correlationOutput <- MicrobiomeDB::selfCorrelation(getCollection(mbioDataset, "16S Genus"), method='spearman', verbose=FALSE)
    correlationDT <- getComputeResult(correlationOutput, "data.table")
    expect_equal(inherits(correlationDT, "data.table"), TRUE)
    expect_equal(all(c('data1', 'data2', 'correlationCoef', 'pValue') %in% names(correlationDT)), TRUE)

    # make sure continuousMetadataOnly flag works so we can do taxa X metadata correlations   
    correlationOutput <- MicrobiomeDB::correlation(genus, method='spearman', verbose=FALSE)
    expect_equal(inherits(correlationOutput, "ComputeResult"), TRUE)

    correlationIGraph <- getComputeResult(correlationOutput, "igraph")
    expect_equal(inherits(correlationIGraph, "igraph"), TRUE)

    correlationBipartiteWidget <- correlationNetwork(correlationOutput)
    expect_equal(inherits(correlationBipartiteWidget, "htmlwidget"), TRUE)
    expect_equal(inherits(correlationBipartiteWidget, "bipartitenetwork"), TRUE)

    selfCorrelationOutput <- MicrobiomeDB::selfCorrelation(genus, method='spearman', verbose=FALSE)
    correlationUnipartiteWidget <- correlationNetwork(selfCorrelationOutput)
    expect_equal(inherits(correlationUnipartiteWidget, "htmlwidget"), TRUE)
    expect_equal(inherits(correlationUnipartiteWidget, "unipartitenetwork"), TRUE)

    # make sure getComputeResultWithMetadata works
    alphaDivOutput <- MicrobiomeDB::alphaDiv(getCollection(mbioDataset, "16S Genus"), method='shannon', verbose=FALSE)
    expect_equal(inherits(alphaDivOutput, "ComputeResult"), TRUE)
    alphaDivDT <- getComputeResultWithMetadata(alphaDivOutput, mbioDataset, metadataVariables = c('country', 'delivery_mode'))
    expect_equal(inherits(alphaDivDT, "data.table"), TRUE)
    expect_equal(all(c('alphaDiversity', 'country', 'delivery_mode') %in% names(alphaDivDT)), TRUE)
})