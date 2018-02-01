


estLocScale <- function(X, type = "wrap", precScale = 1e-12) {
  # Estimates location and scale for every column in X
  #
  # args:
  #   X: data matrix
  #   type: "wrap" or "mcd", the location/scale estimators used
  #   precScale: precision scale used throughout the algorithm
  # Returns: 
  #   loc: the locations of the columns in X
  #   scale: the scales of the columns in X
  
  # Check inputs
  if (!is.data.frame(X) & !is.matrix(X) & !is.vector(X)) {
    stop("The input data must be a vector, matrix or a data frame")
  }
  type <- match(type, c("biwhub","hubhub","wrap","mcd", "rawmcd", "wrapmedmad")) - 1
  if (is.na(type)) {
    stop(paste("Invalid \"type\" argument. Should be \"wrap\" or \"mcd\""))
  }
  
  # Estimate location/scale
  res <- tryCatch( .Call('_cellWise_estLocScale_cpp', as.matrix(X),type,  precScale,
                         PACKAGE = 'cellWise'),
                   "std::range_error" = function(e){
                     conditionMessage( e ) })
  zeroscales <- which(res$scale <= precScale)
  if ( length(zeroscales) > 0) {
    warning(paste(length(zeroscales)," out of ", dim(X)[2], " variables have an estimated scale <= 
\"precScale\" = ", precScale, "."))
  }
  return(list(loc = drop(res$loc), scale = drop(res$scale)))
}