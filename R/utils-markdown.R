# Helper functions for creating the R Markdown files.

#' Create R Markdown file from compute result
#' This function will create an R Markdown file that from a ComputeResult object. The markdown
#' file will describe how to access the data in the ComputeResult object and create a plot. Finally,
#' it will render the R Markdown file to HTML. Templates for these mardown files are stored in the
#' inst/rmarkdown/templates directory.
#' @param result A compute result object
#' @param rmd_file_name The name of the R Markdown file to be created. Default is 'tmp.Rmd'
#' @return A list containing the plots created for the R markdown file.
createRMarkdownFromComputeResult <- function(result, rmd_file_name = "tmp.Rmd"){

  dt <- getComputeResult(result)

  if (result@name == 'betaDiv') {
    # We'll create a pcoa scatterplot
    axis1_name <- result@computedVariableMetadata[[1]]@displayName
    axis2_name <- result@computedVariableMetadata[[2]]@displayName
    p <- ggplot2::ggplot(dt) +
      aes(x=Axis1, y=Axis2) + 
      geom_point() +
      labs(y= axis2_name, x =axis1_name,
            title="Beta diversity by body site",
            caption=paste0("produced on ", Sys.time())) +
      theme_bw()
    
    template <- "inst/rmarkdown/templates/beta_div_mkdn/skeleton/skeleton.Rmd"
    mylist <- list(plot = p)
  } else {
    stop(paste0("This function does not support ComputeResult objects with name=", data@name))
  }

  # Create the R Markdown file
  file.copy(template, rmd_file_name, overwrite = TRUE)
  rmarkdown::render(rmd_file_name)
  
  return(mylist)
  
}


#' Custom HTML template, based on the custom 
#'
#' Loads additional style and template file
#'
#' @param toc should a table of contents be displayed?
#' @param ... additional arguments provided to \@code{html_document}
#' @export
#'
mbiodb_html_format = function(toc = TRUE, ...) {

  # locations of resource files in the package
  pkg_resource = function(...) {
    system.file(..., package = "MicrobiomeDB")
  }

  css    = pkg_resource("rmarkdown/resources/styles.css")

  # call the base html_document function
  rmarkdown::html_document(
    toc = toc,
    toc_float = TRUE,
    fig_width = 6.5,
    fig_height = 4,
    theme = "lumen",
    df_print = "kable",
    number_sections = TRUE,
    ...
  )
}
