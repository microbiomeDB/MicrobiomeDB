<!-- badges: start -->
  [![R-CMD-check](https://github.com/microbiomeDB/MicrobiomeDB/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/microbiomeDB/MicrobiomeDB/actions/workflows/R-CMD-check.yaml)
  [![](https://dcbadge.limes.pink/api/server/DEjNbxgq?style=flat)](https://discord.gg/DEjNbxgq)
  <!-- badges: end -->

# MicrobiomeDB <a href="https://microbiomedb.github.io/MicrobiomeDB/"><img src="man/figures/MicrobiomeDB_hex.png" align="right" height="200" style="float:right; height:200px;"/></a>

A suite of R packages containing all of the data from MicrobiomeDB.org, and tools for analyzing and visualizing the same. It also supports
importing data from the following formats:
 - TreeSummarizedExperiment
 - phyloseq
 - BIOM
 - QIIME2 output files
 - Mothur output files
 - dada2 output files
 - HUMAnN output files

If you have questions or comments of any kind, join us on our [Discord Community Server](https://discord.gg/DEjNbxgq)! We'd love to hear from you.

## Installation

Use the R package [remotes](https://cran.r-project.org/web/packages/remotes/index.html) to install MicrobiomeDB. From the R command prompt:

```{R}
remotes::install_github('microbiomeDB/MicrobiomeDB')
```

If you would like to import data from any format other than TreeSummarizedExperiment, you will also need to install [mia][] from [Bioconductor][] with:

```{R}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("mia")
```

Please note that Linux and MacOS with source-level installations require the 'gsl' system dependency for `mia`. 

On Debian or Ubuntu

```
sudo apt-get install -y libgsl-dev
```

On Fedora, CentOS or RHEL

```
sudo yum install libgsl-devel
```

On macOS (source installations are not common on macOS, so this step is not usually necessary)

```
brew install gsl
```

## Usage
This package is intended to be used to explore the curated datasets from MicrobiomeDB.org, as well as your own datasets. It comes pre-packaged with the same functions used to power the analysis tools from the website. Functions to facilitate easily importing data from phyloseq objects, TreeSummarizedExperiment objects and .biom files (among others) are available.

This package is paired with a dedicated data package called microbiomeData which includes a number of pre-built `MbioDataset` objects representing the curated data from the MicrobiomeDB.org website. You can see their names like:

```R
remotes::install_github('microbiomeDB/microbiomeData')
microbiomeData::getCuratedDatasetNames()
```
This will return a list of names of data objects which were installed with the package. One such dataset is 'DiabImmune', which we'll use in the following example. `MbioDataset` objects contain two things `metadata` and `collections`. Metadata typically include details about samples in the dataset, such as the age of the person they were collected from. Collections can be any group of variables which represent a unified biological concept measured over a consistent or comparable range of values. An example might be relative abundances of various genera. 

```R
getCollectionNames(DiabImmune) # will print the names of collections
myCollection <- getCollection(DiabImmune, '16S Species') # NOTE: you can also use the `format` argument here to get these as phyloseq objects
```

Once you have `myCollection`, you can start using our `microbiomeComputations` package (which was installed for you when you installed this one) to do fun things like:

```R
alphaDivResults <- alphaDiv(myCollection)
correlationResults <- correlation(myCollection)
differentialAbundanceResults <- differentialAbundance(
  myCollection, 
  "breastfed_duration", # see getMetadataVariableNames() and getMetadataVariableSummary()
  groupA = function(x) {x < 300},
  groupB = function(x) {x >= 300},
  method = 'Maaslin2')  
```

This will give you a `ComputeResult` object, with slots for `data` and `statistics` that you can explore. These objects can be difficult to parse, so we've added some functions to help format these results in more usable and exciting ways! They are called `getComputeResult` and `getComputeResultWithMetadata` which will return data.tables (and sometimes igraph objects) which you can use like this:

```R
myCorrelationDT <- getComputeResult(correlationResults)
myCorrelationGraph <- getComputeResult(correlationResults, format = 'igraph')
MicrobiomeDB::correlationNetwork(myCorrelationDT) # will render a network visualization of the results using widgets from our own `corGraph` project

getMetadataVariableNames(DiabImmune) # will print names of metadata variables you can ask for
myAlphaDivDT <- getComputeResultWithMetadata(alphaDivResults, DiabImmune, 'host_body_site')
```

Or you can take these results as a `data.table` object and use them to build plots and things with ggplot2 or any other tool you like. 

Let us know if you build something interesting, encounter any bugs, or just wish something were easier to do. We'd love to hear from you! We consider each new issue filed to be a compliment.

## Contributing
Pull requests are welcome and should be made to the **dev** branch. 

For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0.txt)
