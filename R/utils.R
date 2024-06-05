.require_package <- function(pkg){
    if(!requireNamespace(pkg, quietly = TRUE)){
        stop("'",pkg,"' package not found. Please install the '",pkg,"' package ",
        "to use this feature.", call. = FALSE)
    }
}