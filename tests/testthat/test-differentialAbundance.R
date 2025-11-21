# using DiabImmune for tests here
dataFile1 <- test_path('testdata','DiabImmune/DiabImmune_entity_16SRRNAV4Assay.txt')
metadataFile1 <- test_path('testdata','DiabImmune/DiabImmune_ParticipantRepeatedMeasure.txt')
dataFile2 <- test_path('testdata','DiabImmune/DiabImmune_MetagenomicSequencingAssay.txt')
metadataFile2 <- test_path('testdata','DiabImmune/DiabImmune_Participant.txt')
metadataFile3 <- test_path('testdata','DiabImmune/DiabImmune_Sample.txt')
ontologyFile <- test_path('testdata','DiabImmune/DiabImmune_OntologyMetadata.txt')
mbioDataset <- MbioDataset(list(dataFile1, dataFile2), list(metadataFile2, metadataFile1, metadataFile3), ontologyFile)

genus <- getCollection(mbioDataset, "16S (V4) Genus (Relative taxonomic abundance analysis)")

genusIdCols <- mbioUtils::getIdColumns(genus)
counts <- round(microbiomeComputations::getAbundances(genus, includeIds=FALSE)*1000)
counts <- cbind(genus@data[, genusIdCols, with=FALSE], counts)
genusCounts <- AbsoluteAbundanceData(
    name = 'genusCounts',
    data = counts, 
    sampleMetadata = genus@sampleMetadata, 
    recordIdColumn = genus@recordIdColumn, 
    ancestorIdColumns = genus@ancestorIdColumns
)

test_that("differentialAbundance wrapper works", {    
    # using a variable thats already binary works
    diffAbundOutput <- MicrobiomeDB::differentialAbundance(genus, "delivery_mode", method='Maaslin2', verbose=FALSE)
    expect_equal(inherits(diffAbundOutput, "ComputeResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics, "DifferentialAbundanceResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics@statistics, "data.frame"), TRUE)
    expect_equal(nrow(diffAbundOutput@statistics@statistics) > 0, TRUE)

    # making a variable binary by passing a groupA works
    diffAbundOutput <- MicrobiomeDB::differentialAbundance(genus, "country", groupA = function(x) {x=="Russia"}, method='Maaslin2', verbose=FALSE)
    expect_equal(inherits(diffAbundOutput, "ComputeResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics, "DifferentialAbundanceResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics@statistics, "data.frame"), TRUE)
    expect_equal(nrow(diffAbundOutput@statistics@statistics) > 0, TRUE)

    # making a variable binary by passing groupA and groupB works
    diffAbundOutput <- MicrobiomeDB::differentialAbundance(
        genus, 
        "country", 
        groupA = function(x) {x=="Russia"},
        groupB = function(x) {x %in% c("Finland", "Estonia")},
        method='Maaslin2', 
        verbose=FALSE
    )
    expect_equal(inherits(diffAbundOutput, "ComputeResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics, "DifferentialAbundanceResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics@statistics, "data.frame"), TRUE)
    expect_equal(nrow(diffAbundOutput@statistics@statistics) > 0, TRUE)

    # character vector of values works for categorical vars
    diffAbundOutput <- MicrobiomeDB::differentialAbundance(
        genus, 
        "country", 
        groupA = "Russia",
        groupB = c("Finland", "Estonia"),
        method='Maaslin2', 
        verbose=FALSE
    )
    expect_equal(inherits(diffAbundOutput, "ComputeResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics, "DifferentialAbundanceResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics@statistics, "data.frame"), TRUE)
    expect_equal(nrow(diffAbundOutput@statistics@statistics) > 0, TRUE)

    # character vector of values works for continuous vars errs
    expect_error(
        diffAbundOutput <- MicrobiomeDB::differentialAbundance(
            genus, 
            "breastfed_duration_days", 
            groupA = "100",
            groupB = "300",
            method='Maaslin2', 
            verbose=FALSE
        )
    )

    # excluding some values from both predicates works
    diffAbundOutput <- MicrobiomeDB::differentialAbundance(
        genus, 
        "country", 
        groupA = function(x) {x=="Russia"},
        groupB = function(x) {x %in% c("Finland")},
        method='Maaslin2', 
        verbose=FALSE
    )
    expect_equal(inherits(diffAbundOutput, "ComputeResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics, "DifferentialAbundanceResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics@statistics, "data.frame"), TRUE)
    expect_equal(nrow(diffAbundOutput@statistics@statistics) > 0, TRUE)

    # building predicates on continuous vars works
    diffAbundOutput <- MicrobiomeDB::differentialAbundance(
        genus, 
        "breastfed_duration_days", 
        groupA = function(x) {x<300},
        groupB = function(x) {x>=300},
        method='Maaslin2', 
        verbose=FALSE
    )
    expect_equal(inherits(diffAbundOutput, "ComputeResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics, "DifferentialAbundanceResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics@statistics, "data.frame"), TRUE)
    expect_equal(nrow(diffAbundOutput@statistics@statistics) > 0, TRUE)

    # deseq works as well, so long as you pass it counts
    # thats not easy for users to do currently, but just say they made AbsoluteAbundanceData objects somehow.
    diffAbundOutput <- MicrobiomeDB::differentialAbundance(genusCounts, "delivery_mode", method='DESeq2', verbose=FALSE)
    expect_equal(inherits(diffAbundOutput, "ComputeResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics, "DifferentialAbundanceResult"), TRUE)
    expect_equal(inherits(diffAbundOutput@statistics@statistics, "data.frame"), TRUE)
    expect_equal(nrow(diffAbundOutput@statistics@statistics) > 0, TRUE)

    # passing a variable w 1 value fails
    expect_error(diffAbundOutput <- MicrobiomeDB::differentialAbundance(genus, "host_body_site", method='Maaslin', verbose=FALSE))

    # passing overlapping groupA and groupB fails
    expect_error(
        diffAbundOutput <- MicrobiomeDB::differentialAbundance(
            genus, 
            "breastfed_duration_days", 
            groupA = function(x) {x<300},
            groupB = function(x) {x>=100},
            method='Maaslin2', 
            verbose=FALSE
        )
    )

    # passing only groupA, but where groupA is TRUE for all values fails
    expect_error(
        diffAbundOutput <- MicrobiomeDB::differentialAbundance(
        genus, 
        "breastfed_duration_days", 
        groupA = function(x) {x>=0},
        method='Maaslin2', 
        verbose=FALSE
    )
    )
})

test_that("Maaslin2 wrapper works", {
    # check class of resulting objects, trying random args
    maaslinOutput <- MicrobiomeDB::Maaslin2(
        data = genus, 
        output = tempfile("maaslin"),
        #min_prevalence = 0,
        fixed_effects = 'delivery_mode',
        analysis_method = "LM", # default LM
        normalization = "TSS", # default TSS
        transform = "LOG", # default LOG
        plot_heatmap = F,
        plot_scatter = F)
    # maaslin just dumps a named list, and i dont want to write a test too specific to it in case they change things
    expect_equal(nrow(maaslinOutput$results) > 0, TRUE)

    # make sure we ignore NA values in the fixed_effects variable
    # this should not fail is the main thing
    genus@sampleMetadata@data$delivery_mode[1] <- NA

    maaslinOutput <- MicrobiomeDB::Maaslin2(
        data = genus, 
        output = tempfile("maaslin"),
        #min_prevalence = 0,
        fixed_effects = 'delivery_mode',
        analysis_method = "LM", # default LM
        normalization = "TSS", # default TSS
        transform = "LOG", # default LOG
        plot_heatmap = F,
        plot_scatter = F)
    # maaslin just dumps a named list, and i dont want to write a test too specific to it in case they change things
    expect_equal(nrow(maaslinOutput$results) > 0, TRUE)

    # make sure we ignore empty string values in the fixed_effects variable
    # this should not fail is the main thing
    genus@sampleMetadata@data$delivery_mode[2] <- ''

    maaslinOutput <- MicrobiomeDB::Maaslin2(
        data = genus, 
        output = tempfile("maaslin"),
        #min_prevalence = 0,
        fixed_effects = 'delivery_mode',
        analysis_method = "LM", # default LM
        normalization = "TSS", # default TSS
        transform = "LOG", # default LOG
        plot_heatmap = F,
        plot_scatter = F)
    # maaslin just dumps a named list, and i dont want to write a test too specific to it in case they change things
    expect_equal(nrow(maaslinOutput$results) > 0, TRUE)
})

test_that("DESeq2 wrapper works", {
    # check class of resulting objects, trying random args
    dds <- MicrobiomeDB::DESeqDataSetFromCollection(data = genusCounts,
                                            design = as.formula(paste0("~delivery_mode")),
                                            tidy = FALSE)
    expect_equal(inherits(dds, "DESeqDataSet"), TRUE)

    # make sure you can actually run deseq on resulting object
    # Estimate size factors before running deseq to avoid errors about 0 counts
    geoMeans = apply(DESeq2::counts(dds), 1, function(x){exp(sum(log(x[x > 0]), na.rm=T) / length(x))})
    dds <- DESeq2::estimateSizeFactors(dds, geoMeans = geoMeans)

    # Run DESeq
    deseq_output <- DESeq2::DESeq(dds)
    deseq_results <- DESeq2::results(deseq_output)

    expect_equal(inherits(deseq_results, "DESeqResults"), TRUE)
    expect_equal(nrow(deseq_results) > 0, TRUE)
})
