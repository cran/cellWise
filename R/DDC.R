
DetectDeviatingCells <- function(X, DDCpars = list()){ 
  # X is the input data, and must be a matrix or a data frame.
  #   It must always be provided.
  #
  # DDCpars contains all the options chosen by the user, but has
  # a default:
  #
  if (is.null(DDCpars)) {
    DDCpars <- list()
  }
  if (!is.list(DDCpars)) {
    stop("DDCpars must be a list")
  }
  
  # general parameters
  if (!"precScale" %in% names(DDCpars)) {
    DDCpars$precScale <- 1e-12
  }
  if (!"silent" %in% names(DDCpars)) {
    DDCpars$silent <- FALSE
  }
  if (!"returnBigXimp" %in% names(DDCpars)) {
    DDCpars$returnBigXimp <- FALSE
  }
  
  # parameters for checkDataSet
  if (!"fracNA" %in% names(DDCpars)) {
    DDCpars$fracNA <- 0.5
  }
  if (!"numDiscrete" %in% names(DDCpars)) {
    DDCpars$numDiscrete <- 3
  }
  
  # Parameters for DDC core
  if (!"tolProb" %in% names(DDCpars)) {
    DDCpars$tolProb <- 0.99
  }
  if (!"tolProbCell" %in% names(DDCpars)) {
    DDCpars$tolProbCell <- DDCpars$tolProb
  }
  if (!"tolProbRow" %in% names(DDCpars)) {
    DDCpars$tolProbRow <- DDCpars$tolProb
  }
  if (!"tolProbReg" %in% names(DDCpars)) {
    DDCpars$tolProbReg <- DDCpars$tolProb
  }
  if (!"tolProbCorr" %in% names(DDCpars)) {
    DDCpars$tolProbCorr <- DDCpars$tolProb
  }
  if (!"corrlim" %in% names(DDCpars)) {
    DDCpars$corrlim <- 0.5
  }
  if (!"combinRule" %in% names(DDCpars)) {
    DDCpars$combinRule <- "wmean"
  }
  if (!"includeSelf" %in% names(DDCpars)) {
    DDCpars$includeSelf <- 1
  }
  if (!"rowdetect" %in% names(DDCpars)) {
    DDCpars$rowdetect <- 1
  }
  if (!"numiter" %in% names(DDCpars)) {
    DDCpars$numiter <- 1
  }
  
  # parameters for correlation estimation if fastDDC = TRUE
  if (!"fastDDC" %in% names(DDCpars)) {
    DDCpars$fastDDC <- 1
  }
  if (!"absCorr" %in% names(DDCpars)) {
    DDCpars$absCorr <- 1
  }
  if (!"qdim" %in% names(DDCpars)) {
    DDCpars$qdim <- 20
  }
  if (!"transFun" %in% names(DDCpars)) {
    DDCpars$transFun <- "wrap"
  }
  if (!"bruteForce" %in% names(DDCpars)) {
    DDCpars$bruteForce <- 0
  }
  if (!"treetype" %in% names(DDCpars)) {
    DDCpars$treetype <- "kd"
  }
  if (!"searchtype" %in% names(DDCpars)) {
    DDCpars$searchtype <- "standard"
  }
  if (!"radius" %in% names(DDCpars)) {
    DDCpars$radius <- 0
  }
  if (!"eps" %in% names(DDCpars)) {
    DDCpars$eps <- 0
  }
  if (!"k" %in% names(DDCpars)) {
    DDCpars$k <- 100
  }
  
  # Convert option parameters to integers
  DDCpars$transFun <- pmatch(DDCpars$transFun, 
                             c("Huber", "wrap", "rank"),
                             nomatch = 2)
  DDCpars$treetype <- as.integer(DDCpars$treetype == "bd")
  DDCpars$searchtype <- pmatch(DDCpars$searchtype,
                               c("standard", "priority", "radius"),
                               nomatch = 1)
  DDCpars$combinRule <- pmatch(DDCpars$combinRule,
                               c("wmean","wmedian","mean","median"),
                               nomatch = 1)
  
  # Retrieve parameters from the list:
  #
  fracNA        <- DDCpars$fracNA
  numDiscrete   <- DDCpars$numDiscrete
  precScale     <- DDCpars$precScale
  returnBigXimp <- DDCpars$returnBigXimp
  silent        <- DDCpars$silent
  # Check the data set and set aside columns and rows that do
  # not satisfy the conditions:
  out <- checkDataSet(X, fracNA, numDiscrete, precScale, silent)
  
  checkedData <- out$remX + 0 # +0 puts it in new memory, so cpp doesn't change out$remX
  
  # Carry out the actual DetectDeviatingCells algorithm on
  # the remaining dataset out1$remX :
  res <- tryCatch( .Call('_cellWise_DDC_cpp', checkedData, DDCpars$tolProbCell,
                         DDCpars$tolProbRow, DDCpars$tolProbReg, DDCpars$tolProbCorr,
                         DDCpars$corrlim, DDCpars$combinRule, DDCpars$rowdetect, DDCpars$includeSelf,
                         DDCpars$fastDDC, DDCpars$absCorr, DDCpars$qdim, DDCpars$transFun,
                         DDCpars$treetype, DDCpars$searchtype, DDCpars$radius, DDCpars$eps,
                         DDCpars$bruteForce, DDCpars$k, DDCpars$numiter,
                         DDCpars$precScale,
                         PACKAGE = 'cellWise'),
                   "std::range_error" = function(e){
                     conditionMessage( e ) })
  
  returned.result <- list(Ximp = res$Ximp,k = res$k,ngbrs = res$ngbrs, robcors = res$robcors,
                          robslopes = res$robslopes, Xest = res$Xest, stdResid = res$stdResid,
                          indcells = drop(res$indcells), Ti = res$Ti, indrows = drop(res$indrows),
                          indall = res$indall, indNAs = res$indNAs, Z = res$Z)
  
  if (returnBigXimp) {
    Ximp <- data.matrix(X, rownames.force = TRUE)
    Ximp[out$rowInAnalysis, out$colInAnalysis] <- out$Ximp
    out$Ximp <- Ximp
  }
  return(c(out, returned.result))
}




