test_that("we can get compute results in different formats", {
    mbioDataset <- MicrobiomeDB::DiabImmune
    genus <- getCollection(mbioDataset, "16S Genus", continuousMetadataOnly = TRUE)

    # make sure metadata dont contain IRIs
    expect_equal(all(grepl('[',names(genus@sampleMetadata@data),fixed=T)), FALSE)
    # and that that means diff abund works now
    comparatorVariable <- microbiomeComputations::Comparator(
                        variable = veupathUtils::VariableMetadata(
                            variableSpec = vuepathUtils::VariableSpec(
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
    diffAbundOutput <- MicrobiomeDB::differentialAbundance(getCollection(mbioDataset, "16S Genus"), comparatorVariable, method='Maaslin', verbose=FALSE)
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

    # make sure getComputeResultWithMetadata works
    alphaDivOutput <- MicrobiomeDB::alphaDiv(getCollection(mbioDataset, "16S Genus"), method='shannon', verbose=FALSE)
    expect_equal(inherits(alphaDivOutput, "ComputeResult"), TRUE)
    alphaDivDT <- getComputeResultWithMetadata(alphaDivOutput, mbioDataset, metadataVariables = c('country', 'delivery_mode'))
    expect_equal(inherits(alphaDivDT, "data.table"), TRUE)
    expect_equal(all(c('alphaDiversity', 'country', 'delivery_mode') %in% names(alphaDivDT)), TRUE)
})