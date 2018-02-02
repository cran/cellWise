
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
  DDCpars$transFun <- pmatch(DDCpars$transFun, 
                             c("Huber", "wrap", "rank"))
  DDCpars$treetype <- as.integer(DDCpars$treetype == "bd")
  DDCpars$searchtype <- pmatch(DDCpars$searchtype,
                               c("standard", "priority", "radius"))
  DDCpars$combinRule <- pmatch(DDCpars$combinRule, c("wmean","wmedian","mean","median"))
  # 
  # Meaning of these options:
  # fracNA     : only consider columns and rows with fewer NAs (missing
  #              values) than this fraction (percentage).
  # numDiscrete: a column that takes on numDiscrete or fewer values will
  #              be considered discrete and not used in the analysis.
  # precScale  : only consider columns whose scale is > precScale.
  #              Here scale is measured by the median absolute deviation.  
  # tolProb    : tolerance probability, with default 0.99, which
  #              determines the cutoff values for flagging outliers.
  #              Used in all the steps of the algorithm.
  # corrlim    : when tring to estimate z_ij from other variables h, we 
  #              will only use variables h with abs(corr(j,h)) >= corrlim.
  #              Variables j without any correlated variables h satisfying 
  #              this are considered standalone, and treated on their own.  
  # combinRule : the operation to combine estimates of z_ij coming from
  #              other variables h: can be wmean, wmedian, mean, median.
  # includeSelf: whether or not the combination rule will include the
  #              variable j itself.
  # rowdetect   : whether the rule for flagging rows is to be applied.
  # returnBigXimp : if TRUE, the imputed data matrix Ximp in the output
  #                 will include the rows and columns that were not
  #                 part of the analysis (and can still contain NAs).
  #
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
  res <- tryCatch( .Call('_cellWise_DDC_cpp', checkedData, DDCpars$tolProb, DDCpars$corrlim,
                         DDCpars$combinRule, DDCpars$rowdetect, DDCpars$includeSelf,
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




