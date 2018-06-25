

wrap <- function(X, locX, scaleX, precScale = 1e-12) {
  # Transforms multivariate data X using wrapping function
  # with b = 1.5 and c = 4 and loc&scale in
  # the input arguments locX and scaleX
  #
  # args:
  #   X: data matrix
  #   locX: vector with location estimates for every column of X
  #   scaleX: vector with scale estimates for every column of X
  #   precScale: precision scale used throughout the algorithm
  # Returns: 
  #   Xw: a matrix with wrapped data
  #   colInWrap: vector of indices of wrapped columns
  #
  
  # Check inputs
  if (!is.data.frame(X) & !is.matrix(X) & !is.vector(X)) {
    stop("The input data must be a vector, matrix or a data frame")
  }
  X <- as.matrix(X)
  if (!length(locX) == dim(X)[2] || !length(scaleX) == dim(X)[2]) {
    stop("The arguments \"locX\" and \"scaleX\" should both be of length dim(X)[2]")
  }
  colInWrap <- which(scaleX > precScale)
  if (length(colInWrap) == 0) {
    stop(paste("No columns with scale > precScale were found."))
  }

  if (precScale < 1e-12) {
    precScale <- 1e-12
    warning("The \"precScale\" argument was changed to 1e-12.")
  }
  
  # Execute wrapping
  res <- tryCatch( .Call('_cellWise_Wrap_cpp', X[, colInWrap], locX, scaleX, precScale,
                         PACKAGE = 'cellWise'),
                   "std::range_error" = function(e){
                     conditionMessage( e ) })
  
  if (length(colInWrap) < dim(X)[2]) {
    warning(paste(dim(X)[2] - length(colInWrap),
                  " variable(s) had a scale <= precScale and were left out.\n
                    The indices of the remaining columns are in $colInWrap."))
  }
  
  return(list(Xw = res$Xw, colInWrap = colInWrap))
}