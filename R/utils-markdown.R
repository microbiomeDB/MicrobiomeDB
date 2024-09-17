# Helper functions for creating the R Markdown files.


createRMarkdownFromComputeResult <- function(data, rmd_file_name = "tmp.Rmd"){

  df <- getComputeResult(data)
  axis1_name <- data@computedVariableMetadata[[1]]@displayName
  axis2_name <- data@computedVariableMetadata[[2]]@displayName
  p <- ggplot2::ggplot(df) +
    aes(x=Axis1, y=Axis2) + 
    geom_point() +
    labs(y= axis2_name, x =axis1_name,
          title="Beta diversity by body site",
          caption=paste0("produced on ", Sys.time())) +
    theme_bw()
  
  template <- "inst/rmarkdown/templates/beta_div_mkdn/skeleton/skeleton.Rmd"
  file.copy(template, rmd_file_name, overwrite = TRUE)
  rmarkdown::render(rmd_file_name)
  
  mylist <- list(plot = p, df = df)
  return(mylist)
  
}
