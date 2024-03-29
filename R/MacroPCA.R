
MacroPCA <- function(X, k = 0, MacroPCApars = NULL) {
  #
  # This function performs the MacroPCA algorithm.
  #
  # Its inputs are:
  #
  # X : the input data, which must be a matrix or a data frame.
  #     It must always be provided.
  # k : the desired number of principal components. 
  #     If k = 0 or k = NULL, the algorithm will compute the
  #     percentage of explained variability for k upto kmax and
  #     show a scree plot, and suggest to choose a value of
  #     k such that the cumulative percentage of explained
  #     variability is >= 80% .
  # MacroPCApars : a list with all the options chosen by the user.
  #     Its default is
  #     MacroPCApars <- list(DDCpars=NULL,kmax=10,alpha=0.50,
  #                          h=NULL,scale=TRUE,maxdir=250,
  #                          distprob=0.99,silent = TRUE, 
  #                          maxiter=20, tol=0.005, center = NULL)
  #
  # Meaning of these options:
  # DDCpars   : list with parameters for the first step of the 
  #             MacroPCA algorithm (for the complete list see the 
  #             function DDC). Default is NULL.
  # kmax      : maximal number of principal components to compute.
  #             Default is kmax=10. 
  #             If k is provided kmax does not need to be specified,
  #             unless k is larger than 10 in which case you need
  #             to set kmax high enough.
  # alpha     : this is the coverage, i.e. the fraction of rows the
  #             algorithm should give full weight. Alpha should be
  #             between 0.50 and 1, the default is 0.50.
  # scale     : a value indicating whether and how the original 
  #             variables should be scaled. If scale=FALSE (default) 
  #             or scale=NULL no scaling is performed (and a vector 
  #             of 1s is returned in the $scaleX slot).
  #             If scale=TRUE the data are scaled to have a robust
  #             scale of 1.
  #             Alternatively scale can be a function like mad,
  #             or a vector of length equal to the number of columns 
  #             of x. 
  #             The resulting scale estimates are returned in the
  #             $scaleX slot of the MacroPCA output.
  # maxdir    : maximal number of random directions to use for
  #             computing the outlyingness of the data points. 
  #             Default is maxdir=250. If the number n of observations
  #             is small all n*(n-1)/2 pairs of observations are used.
  # distprob  : quantile determining the cutoff values for
  #             orthogonal and score distances. Default is 0.99.
  # silent     : whether to print intermediate results. 
  #             Default is TRUE 
  # maxiter   : maximum number of iterations. Default is 20.
  # tol       : tolerance for iterations. Default is 0.005.
  # bigOutput : whether to compute and return NAimp, Fullimp and Cellimp
  # center    : if NULL, MacroPCA will compute the center. If a vector 
  #             with d components, this center will be used.
  #
  # The outputs are:
  #
  # MacroPCApars: the options used in the call.
  # scaleX      : the scales of the columns of X.
  # k           : the number of principal components.  
  # loadings    : the columns are the k loading vectors.
  # eigenvalues : the k eigenvalues.
  # center      : vector with the fitted center.
  # alpha       : alpha from the input.
  # h           : h (computed from alpha).
  # It          : number of iteration steps.
  # diff        : convergence criterion.
  # X.NAimp     : data with all NA's imputed by MacroPCA.
  # scores      : scores of X.NAimp 
  # OD          : orthogonal distances of the rows of X.NAimp 
  # cutoffOD    : cutoff value for the OD.
  # highOD      : row numbers of observations with OD > cutoffOD.
  # SD          : score distances of the rows of X.NAimp 
  # cutoffSD    : cutoff value for the SD.
  # highSD      : row numbers of observations with SD > cutoffSD.
  # residScale  : scale of the residuals.
  # stdResid    : standardized residuals. Note that these are NA
  #               for all missing values of the data X.
  # indcells    : indices of cellwise outliers.
  # NAimp       : various results for the NA-imputed data.
  # Cellimp     : various results for the cell-imputed data.
  # Fullimp     : various result for the fully imputed data.
  # DDC         : results of the first step of MacroPCA. These are
  #               needed to run MacroPCApredict on new data.
  
  # The random seed is retained when leaving the function
  if (exists(".Random.seed",envir = .GlobalEnv,inherits = FALSE)) {
    seed.keep <- get(".Random.seed", envir = .GlobalEnv, 
                     inherits = FALSE)
    on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
  }
  
  # Check input parameters:
  if (!is.data.frame(X) & !is.matrix(X)) {
    stop("The input data must be a matrix or a data frame.")
  }
  d <- ncol(X)
  if (d < 2) stop("The input data must have at least 2 columns.")
  
  if (is.null(MacroPCApars)) {
    MacroPCApars <- list()
  }
  if (!is.list(MacroPCApars)) {
    stop("MacroPCApars must be a list")
  }
  
  # Parameters for checkDataset
  if (!"DDCpars" %in% names(MacroPCApars)) {
    MacroPCApars$DDCpars <- list()
  }
  if (!"fracNA" %in% names(MacroPCApars$DDCpars)) {
    MacroPCApars$DDCpars$fracNA <- 0.5
  }
  if (!"numDiscrete" %in% names(MacroPCApars$DDCpars)) {
    MacroPCApars$DDCpars$numDiscrete <- 3
  }
  if (!"cleanNAfirst" %in% names(MacroPCApars$DDCpars)) {
    MacroPCApars$DDCpars$cleanNAfirst <- "automatic"
  }
  if (!"silent" %in% names(MacroPCApars$DDCpars)) {
    MacroPCApars$DDCpars$silent <- if (is.null(MacroPCApars$silent)) {
      FALSE
    } else {
      MacroPCApars$silent
    }
  }
  # We want the ability to set MacroPCApars$DDCpars$center
  # separately from MacroPCApars$center :
  if (!("center" %in% names(MacroPCApars$DDCpars))) {
    DDCcenter <- NA
  } else {
    DDCcenter = MacroPCApars$DDCpars$center
    if(is.na(DDCcenter[1])){
      if (length(DDCcenter) > 1) {
        stop(paste0("MacroPCApars$DDCpars$center should be NA,\n",
                    "or a vector of length ncol(X) without NAs"))
      }
    } else {
      DDCcenter <- as.vector(DDCcenter)
      ddc = length(DDCcenter)
      print(dim(X))
      if(ddc != d) stop(paste0(
        "length(MacroPCApars$DDCpars$center) = ",ddc,
        " should equal ncol(X) = ",d))
      if (any(is.na(DDCcenter))) stop(
        "MacroPCApars$DDCpars$center should not contain NAs")
    }
  } 
  MacroPCApars$DDCpars$center <- DDCcenter
  
  # parameters for MacroPCA
  if (!"kmax" %in% names(MacroPCApars)) {
    MacroPCApars$kmax <- 10
  }
  if (!"alpha" %in% names(MacroPCApars)) {
    MacroPCApars$alpha <- 0.50
  }
  if (!"scale" %in% names(MacroPCApars)) {
    MacroPCApars$scale <- TRUE
  }
  if (!"maxdir" %in% names(MacroPCApars)) {
    MacroPCApars$maxdir <- 250
  }
  if (!"distprob" %in% names(MacroPCApars)) {
    MacroPCApars$distprob <- 0.99
  }
  if (!"silent" %in% names(MacroPCApars)) {
    MacroPCApars$silent <- FALSE
  }
  if (!"maxiter" %in% names(MacroPCApars)) {
    MacroPCApars$maxiter <- 20
  }
  if (!"tol" %in% names(MacroPCApars)) {
    MacroPCApars$tol <- 0.005
  }
  if (!"bigOutput" %in% names(MacroPCApars)) {
    MacroPCApars$bigOutput <- TRUE
  }
  if (!"center" %in% names(MacroPCApars)) {
    MacroPCApars$center <- NA
  }
  
  # Put options in convenient order
  MacroPCApars <- MacroPCApars[c("DDCpars", "kmax", "alpha", 
                                 "scale", "maxdir", "distprob",
                                 "silent", "maxiter", "tol",
                                 "bigOutput", "center")]
  
  # Input arguments
  DDCpars   <- MacroPCApars$DDCpars
  kmax      <- MacroPCApars$kmax
  alpha     <- MacroPCApars$alpha
  maxdir    <- MacroPCApars$maxdir
  distprob  <- MacroPCApars$distprob
  silent    <- MacroPCApars$silent
  maxiter   <- MacroPCApars$maxiter
  tol       <- MacroPCApars$tol
  scale     <- MacroPCApars$scale
  bigOutput <- MacroPCApars$bigOutput
  center    <- MacroPCApars$center
  
  # check fixedCenter
  if(is.na(center[1])){
    fixedCenter <- FALSE
    if (length(center) > 1) {
      stop("MacroPCApars$center should be a vector of length ncol(X) without NAs")
    }
  } else {
    fixedCenter <- TRUE
    givenCenter <- as.vector(center)
    dc = length(givenCenter)
    if(dc != d) stop(paste0(
      "length(MacroPCApars$center) = ",dc," should equal ncol(X) = ",d))
    if (any(is.na(center))) stop("MacroPCApars$center should not contain NAs")
  }
  
  # check dataset (some rows/columns may leave the analysis).
  checkOut <- cellWise::checkDataSet(X, DDCpars$fracNA,
                                     DDCpars$numDiscrete,
                                     DDCpars$precScale, 
                                     DDCpars$silent,
                                     DDCpars$cleanNAfirst)
  X <- as.matrix(checkOut$remX)
  n <- nrow(X)
  d <- ncol(X)
  if(fixedCenter == TRUE){
    givenCenter <- givenCenter[checkOut$colInAnalysis]
    X <- scale(X, center = givenCenter, scale = FALSE)
    # drop attribute center:
    attr(X, "scaled:center") <- NULL
  }
  if(!is.na(DDCcenter[1])) DDCcenter <- DDCcenter[checkOut$colInAnalysis] 
  
  ## Verify inputs
  
  if (is.null(k)) { k <- 0 }
  
  # Check kmax, k, alpha, h
  if (kmax < 1) stop("kmax should be at least 1.")
  kmax <- min(kmax, d) # PR
  k <- floor(k)
  if (k < 0) {
    k <- 0
  } else if (k > kmax) {
    cat(paste("\nThe number of principal components k = ", 
              k, " is larger than kmax = ", kmax,
              "\n so k is set to ", kmax, ".\n", sep = ""))
    k <- kmax
  }
  
  if (alpha < 0.5 | alpha > 1) 
    stop("Alpha is out of range: should be between 1/2 and 1")
  if (k == 0) 
    h <- h.alpha.n(alpha, n, kmax)
  else h <- h.alpha.n(alpha, n, k)
  # cat(paste0("h = ",h,"\n"))
  
  
  ##########################
  ##  Detect deviating cells
  ##########################
  
  DDCpars$coreOnly <- TRUE 
  ## now only execute DDCcore, since CheckDataSet already done
  if(!is.na(DDCcenter[1])) DDCpars$center <- DDCcenter
  resultDDC <- cellWise::DDC(X, DDCpars)
  DDCpars   <- resultDDC$DDCpars # does not contain the center
  if(!is.na(DDCcenter[1])) DDCpars$center <- DDCcenter
  MacroPCApars$DDCpars <- DDCpars
  
  Ti <- resultDDC$Ti
  indrows1 <- which(abs(Ti) > sqrt(qchisq(0.99, 1)))
  indrows2 <- intersect(order(abs(Ti),decreasing = T)[seq_len(n - h)],
                        indrows1)
  indcells <- resultDDC$indcells    
  XO       <- X   # original matrix, with its NA's
  Xnai     <- X   # initialize NA-imputed data matrix
  Xci      <- X   # initialize cell-imputed data matrix
  indNA    <- which(is.na(X))            
  # positions of missing values
  indimp   <- unique(c(indcells, indNA))
  # positions where values can be imputed
  Xci[indimp] <- resultDDC$Ximp[indimp] 
  # Imputes all missings and cellwise outliers
  # This is not the same as resultDDC$Ximp, 
  # which imputes indrows completely.
  
  
  #########################
  # Step 1: Standardization
  #########################
  
  # Compute scales of the variables if requested
  if (is.logical(scale)) {
    if (!scale) scaleX <- rep(1, d) # vector("numeric", d) + 1
    if (scale) {
      locScale <- cellWise::estLocScale(X, nLocScale = DDCpars$nLocScale,
                                        type = "1stepM", 
                                        center = !fixedCenter)
      scaleX <- locScale$scale
    }
  } else if (is.function(scale)) {
    scaleX <- apply(X, 2, FUN = scale, na.rm = TRUE)
  } else if (is.vector(scale)) scaleX <- scale  
  rm(X) # to save space
  
  # Scale the variables
  XO   <- sweep(XO, 2, scaleX, "/")
  Xnai <- sweep(Xnai, 2, scaleX, "/")  
  Xci  <- sweep(Xci, 2, scaleX, "/")  
  # Note that scaleX is 1 when scale=F
  
  ############################
  # Step 2: Projection pursuit
  ############################
  
  # Compute initial cell-imputed matrix
  XOimp <- Xci  # store imputed rowoutliers    
  Xind  <- Xci; Xind[indimp] <- NA; Xind <- is.na(Xind)
  
  # impute original NA values in unimputed rowoutliers
  Xnai[indNA] <- XOimp[indNA] 
  
  # classical PCA to determine the rank of the data matrix
  # (this step is the same as in ROBPCA)
  XciSVD <- truncPC(Xci, center = !fixedCenter)
  
  if (XciSVD$rank == 0) stop("All data points collapse!")
  
  center <- rep(0, d)  
  rot    <- diag(d)
  
  rowcellind   <- tabulate(((indcells - 1) %% n) + 1, n)
  rowcellind[indrows2] <- d # pretends all cells flagged, so
  # the rows flagged as outlying by DDC will not be imputed.
  hrowsLowcell <- order(rowcellind)[seq_len(h)]
  indComb <- unique(hrowsLowcell,
                    which(rowcellind == rowcellind[hrowsLowcell[h]]))
  # these rows have the fewest cellwise outliers (includes ties)
  indComb <- setdiff(seq_len(n), indComb) # take complement
  
  # Use unimputed rowwise outliers for least outlying points 
  # (original NA values are imputed)
  
  if (length(indComb) > 0) Xci[indComb, ] <- Xnai[indComb, ]
  
  # Find least outlying points by projections:
  
  if (!fixedCenter) {
    alldir <- choose(n, 2)
    ndir <- min(maxdir, alldir) # too low for fixedCenter
    all  <- (ndir == alldir) # <==> (maxdir >= alldir)
    B <- makeDirections(Xci, ndir, all = all,
                         fixedCenter = FALSE)
    # Calculates directions through two data points.
  } else {
    Xnorm <- apply(Xci, 1, vecnorm)
    nzind  <- which(Xnorm > 1e-12) # to avoid dividing by zero
    ndir <- maxdir # we will always use maxdir directions
    B <- makeDirections(Xci[nzind, ], ndir, all = FALSE,
                         fixedCenter = TRUE)
    # Calculates ndir directions between the origin and one 
    # datapoint, and if there are not enough of them it also 
    # uses averages of some of these unit vectors.
  }
  Bnorm <- apply(B, 1, vecnorm)
  nonz <- which(Bnorm > 1e-12) # to avoid dividing by zero
  Bnormr <- Bnorm[nonz]
  m <- length(nonz)
  B <- B[nonz, ]
  A <- diag(1 / Bnormr) %*% B
  Y <- Xci %*% t(A) # projected points in columns
  Z <- matrix(0, n, m)
  
  umcds <- cellWise::estLocScale(Y,nLocScale = DDCpars$nLocScale,
                                 type = "mcd", alpha = alpha,
                                 center = !fixedCenter,
                                 silent = TRUE)
  
  zeroscales <- which(umcds$scale < 1e-12)
  for (i in seq_len(length(zeroscales))) {
    umcdweights <- unimcd(Y[, zeroscales[i]], alpha = alpha,
                          center = !fixedCenter)$weights
    if (robustbase::rankMM(Xci[umcdweights == 1, ]) == 1) {
      stop("At least ", sum(umcdweights),
           " observations are identical.")
    }
  }
  
  # There should be no more zero scales since we didn't stop.
  Z <- abs(scale(Y, umcds$loc, umcds$scale))
  
  H0 <- order(apply(Z, 1, max))
  H0 <- setdiff(H0,indrows2)
  H0 <- H0[seq_len(h)]
  
  
  ############################
  # Step 3: Subspace dimension
  ############################  
  
  # Use imputed rowwise outliers for reconstruction
  if (length(indComb) > 0) Xci[indComb, ] <- XOimp[indComb, ]
  
  # determine k
  Xcih <- Xci[H0, ]
  Xcih.SVD <- truncPC(Xcih, ncomp = kmax, center = !fixedCenter)
  kmax <- min(Xcih.SVD$rank, kmax)
  
  if (!silent) {
    printvals = Xcih.SVD$eigenvalues
    if(length(printvals) < d)
      printvals = c(Xcih.SVD$eigenvalues, 0)
    if(length(printvals) < d){
      cat("\nInitial eigenvalues:\n", printvals, " ... \n")
    } else {
      cat("\nInitial eigenvalues:\n", printvals, "\n")
    }
  }
  
  if (k == 0) {
    # To help the user determine the number of PC's, we
    # compute cumulative variability and scree plot and stop.
    # This necessitates computing _all_ eigenvalues:
    Xcih.SVD <- truncPC(Xcih, ncomp = min(n, d),
                        center = !fixedCenter)
    ratios = Xcih.SVD$eigenvalues / Xcih.SVD$eigenvalues[1]
    test <- which(ratios <= 0.001)
    k <- if (length(test) != 0) 
      min(min(Xcih.SVD$rank, test[1]), kmax)
    else min(Xcih.SVD$rank, kmax)
    cumulative <- cumsum(Xcih.SVD$eigenvalues[seq_len(k)]) /
      sum(Xcih.SVD$eigenvalues)
    # We used _all_ eigenvalues in the denominator.
    barplot(Xcih.SVD$eigenvalues[seq_len(k)], main = "scree plot",
            ylab = "eigenvalues", names.arg = seq_len(k))
    
    roundcumul = round(100 * cumulative, 1)
    names(roundcumul) = paste("PC", seq_len(length(cumulative)),sep = "")
    cat(paste(c("\nThe cumulative percentage of explained variability",
                "\nby the first", k, "components is:\n")))
    print(noquote(format(roundcumul, nsmall=1))) 
    
    if (cumulative[k] > 0.8) {
      kselect <- which(cumulative >= 0.8)[1]
      cat(paste("\nBased on explained variability >= 80% one would",
                " select k = ", kselect, ".\n", sep = ""))
    } else {
      cat(paste("\nThe selected value of kmax does not allow for an explained ",
                "\nvariance of at least 80%, since kmax components only explain\n",
                roundcumul[k],
                "% of variance. Consider increasing the value of kmax."))
    }
    cat(paste("\nPlease use this information and the scree plot",
              " to select a value of k",
              "\nand rerun MacroPCA with it.\n\n", sep = "")) 
    return(list(eigenvalues = Xcih.SVD$eigenvalues[seq_len(k)],
                cumulativeVar = cumulative))
  }
  
  if (!silent)
    cat("\nXciSVD$rank, Xcih.SVD$rank, k and kmax: ", XciSVD$rank,
        Xcih.SVD$rank, k, kmax, "\n")
  
  
  #######################################
  # Step 4: Iterative subspace estimation
  #######################################    
  
  It     <- 0 # this is the iteration index s in the paper
  diff   <- 0
  k      <- min(k, Xcih.SVD$rank)
  Pr     <- Xcih.SVD$loadings[, seq_len(k)]
  PrPrev <- Pr
  mXci   <- Xcih.SVD$center # center
  if (any(Xind) & maxiter > 0) { 
    # Step3: Iterate
    diff <- 100
    while (It < maxiter & diff > tol) { # Iterate
      It        <- It + 1;
      XciC      <- sweep(Xci, 2, mXci) # centered Xci  
      Tr        <- XciC %*% Pr         # scores matrix     
      Xcihat    <- Tr %*% t(Pr)        # fit to XciC                     
      Xcihat    <- sweep(Xcihat, 2, mXci, "+") # fit to Xci  
      Xci[Xind] <- Xcihat[Xind]         
      # impute missings + cellwise outliers
      
      Xcih <- Xci[H0, ]
      Xcih.SVD <- truncPC(Xcih, ncomp = k,
                          center = !fixedCenter)
      k <- min(k, Xcih.SVD$rank)
      Pr <- Xcih.SVD$loadings[, seq_len(k)]  # loadings matrix
      mXci <- Xcih.SVD$center        # mean vector
      
      diff <- maxAngle(Pr, PrPrev)
      PrPrev <- Pr
    } # Iteration ends
  } 
  
  XciC      <- sweep(Xci, 2, mXci) # centered Xci  
  Tr        <- XciC %*% Pr         # scores matrix
  Xcihat    <- Tr %*% t(Pr)        # fit to XciC                
  Xcihat    <- sweep(Xcihat, 2, mXci, "+")  # fit to Xci 
  Xci[Xind] <- Xcihat[Xind]         
  # impute missings + cellwise
  
  Xnai[indNA] <- Xci[indNA] 
  # update imputations of missings in Xnai
  indrows3    <- setdiff((seq_len(n)), H0) 
  # indrows3 means: all except H0
  
  # initialize Xfi, the fully imputed X
  Xfi <- Xci
  
  # Xci should not impute outlying cells in outlying rows:
  if (length(indrows3) > 0) {
    Xci[indrows3, ] <- Xnai[indrows3, ]
  }
  
  
  #####################
  # Step 5: reweighting
  #####################
  
  if (k < XciSVD$rank) {  
    if (!silent) 
      cat("\nPerformed an extra reweighting step because k =", k,
          "< rank =", XciSVD$rank, ".\n")
    XciRc  <- Xci - matrix(rep(mXci, times = n), 
                           nrow = n, byrow = TRUE)
    Xcihat <- XciRc %*% Pr %*% t(Pr) # fit to XciRc
    Rdiff  <- XciRc - Xcihat # orthogonal residual
    ODh    <- apply(Rdiff, 1, vecnorm) # norm of residual
    umcd   <- unimcd(ODh ^ (2 / 3), alpha) 
    # its robust location and scale
    cutoffODh <- sqrt(qnorm(distprob, umcd$loc, umcd$scale) ^ 3)
    # This is the cutoff on the orthogonal distances of the
    # cell-imputed data.
    
    Hstar    <- (ODh <= cutoffODh) 
    # The set Hstar is denoted $H^*$ in the paper.
    Hstar[indrows2] <- FALSE 
    # removes outlying rows flagged by DDC
    Xcih.SVD <- truncPC(Xfi[Hstar, ], ncomp = k,
                        center = !fixedCenter)
    k <- min(Xcih.SVD$rank, k)
  } else {
    Hstar     <- rep(0, dim(Xci)[1])
    Hstar[H0] <- 1
    Hstar     <- as.logical(Hstar)
  }
  # Hstar is the final set of rowwise inliers.  
  
  indrows4 <- which(!Hstar) # complement
  
  # Construct the final cell-imputed Xci
  Xci <- Xfi
  if (length(indrows4) > 0) {
    # In the outlying rows we only impute the NAs:
    Xci[indrows4, ] <- Xnai[indrows4, ]
  }
  
  # compute new center $\mu^*$ and new loading matrix $P^*$:
  center <- center + Xcih.SVD$center %*% t(rot) # new center,
  # but in case of fixedCenter the new center remains zero.
  rot  <- rot %*% Xcih.SVD$loadings
  n1   <- sum(Hstar)
  h1   <- h.alpha.n(alpha, n1, k)
  Xci1 <- (Xci[Hstar, ] - matrix(rep(Xcih.SVD$center, times = n1), 
                                 nrow = n1, byrow = TRUE)) %*% 
    Xcih.SVD$loadings
  Xci1 <- as.matrix(Xci1[, seq_len(k)])
  rot  <- as.matrix(rot[, seq_len(k)]) # new loading matrix
  mah  <- mahalanobis(Xci1, center = rep(0, ncol(Xci1)), 
                      cov = diag(Xcih.SVD$eigenvalues[seq_len(k)], 
                                 nrow = k))
  
  
  ###############
  # Step6: DetMCD
  ###############
  
  oldobj <- prod(Xcih.SVD$eigenvalues[seq_len(k)])
  niter  <- 100
  for (j in seq_len(niter)) { # This part is from ROBPCA
    Xcih     <- as.matrix(Xci1[order(mah)[seq_len(h1)], ], ncol = k)
    Xcih.SVD <- truncPC(Xcih, ncomp = k,
                        center = !fixedCenter)
    obj    <- prod(Xcih.SVD$eigenvalues)
    Xci1   <- (Xci1 - matrix(rep(Xcih.SVD$center, times = n1), 
                             nrow = n1, byrow = TRUE)) %*% Xcih.SVD$loadings
    center <- center + Xcih.SVD$center %*% t(rot)
    rot    <- rot %*% Xcih.SVD$loadings
    mah    <- mahalanobis(Xci1, center = matrix(0, 1, ncol(Xci1)), 
                          cov = diag(Xcih.SVD$eigenvalues, 
                                     nrow = length(Xcih.SVD$eigenvalues)))
    if (Xcih.SVD$rank == k & abs(oldobj - obj) < 1e-12) 
      break
    oldobj <- obj
    if (Xcih.SVD$rank < k) {
      j <- 1
      k <- Xcih.SVD$rank
    }
  }
  
  if (fixedCenter) {
    Xci2mcd <- covMcd2(Xci1, nsamp = "deterministic", alpha = h1 / n1,
                       center = rep(0, ncol(Xci1)))
  } else {
    Xci2mcd <- rrcov::CovMcd(Xci1, nsamp = "deterministic", alpha = h1 / n1)
  }
  
  eps <- 1e-16
  crit <- if (fixedCenter) {prod(eigen(Xci2mcd$cov, only.values = TRUE)$values)} else {Xci2mcd@crit}
  if (crit < obj + eps) {
    if (fixedCenter) {
      Xci2cov    <- Xci2mcd$cov
      Xci2center <- Xci2mcd$center
      
    } else {
      Xci2cov    <- rrcov::getCov(Xci2mcd)
      Xci2center <- rrcov::getCenter(Xci2mcd)
    }
    if (!silent) 
      cat("\nThe final step used eigenvectors of MCD scatter.\n")
  } else {
    consistencyfactor <- median(mah) / qchisq(0.5, k)
    mah     <- mah / consistencyfactor
    weights <- ifelse(mah <= qchisq(0.975, k), TRUE, 
                      FALSE)
    wcov    <- cov.wt(x = Xci1, wt = weights, center = !fixedCenter,
                      method = "ML")
    if (fixedCenter) {wcov$center <- rep(0, ncol(Xci1))}
    Xci2center <- wcov$center
    Xci2cov    <- wcov$cov
    if (!silent) 
      cat("\nThe final step used eigenvectors of reweighted covariance.\n")
  }
  ee     <- eigen(Xci2cov)
  P6     <- ee$vectors
  center <- as.vector(center + Xci2center %*% t(rot))
  eigenvalues <- ee$values
  loadings    <- rot %*% P6
  # flip signs of final loadings to make result more unique:
  flipcolumn = function(x){ 
    if (x[which.max(abs(x))] < 0) {-x} else {x} }
  loadings = apply(loadings, 2L, FUN = flipcolumn)
  
  dimnames(loadings) <- list(colnames(Xnai), paste(
    "PC", seq_len(ncol(loadings)), sep = ""))
  if(!fixedCenter){
    MacroPCApars["center"] <- NULL # to match pars of old MacroPCA
  }
  # Start the output list: res
  res <- list(MacroPCApars = MacroPCApars, remX = checkOut$remX,
              DDC = c(checkOut, resultDDC), scaleX = scaleX, k = k,
              loadings = loadings, eigenvalues = eigenvalues,
              center = center, alpha = alpha, h = h,
              It = It, diff = diff) 
  # All of these elements are common to Xci, Xnai, Xfi
  
  
  #####################################
  # Step 7: Scores and predicted values
  #####################################
  
  # cell-imputed matrix Xci : scores and distances
  
  scoresci <- (Xci - matrix(rep(center, times = n), nrow = n, 
                            byrow = TRUE)) %*% loadings
  Cellimp  <- list(scoresci = scoresci)
  out      <- pca.distancesNew(res, Xci, scoresci,
                               XciSVD$rank, distprob)
  names(out)[1] <- "ODci"
  names(out)[3] <- "SDci"
  cutoffOD      <- out$cutoffOD
  out$highODci  <- which(out$ODci > cutoffOD)
  cutoffSD      <- out$cutoffSD
  out$highSDci  <- which(out$SD > cutoffSD)
  Cellimp       <- c(Cellimp, out); rm(out)
  
  # NA-imputed data Xnai : scores and distances
  res$X.NAimp <- Xnai
  scoresnai   <- (Xnai - matrix(rep(center, times = n), nrow = n, 
                                byrow = TRUE)) %*% loadings
  res$scores  <- scoresnai
  
  NAimp          <- list(scoresnai = scoresnai)
  out            <- pca.distancesNew(res, Xnai, scoresnai,
                                     XciSVD$rank, distprob)
  res$OD         <- out$OD
  out$cutoffOD   <- cutoffOD # cutoff based on the cell-imputed
  res$cutoffOD   <- out$cutoffOD
  out$highODnai  <- which(out$OD > cutoffOD)
  res$SD         <- out$SD
  out$cutoffSD   <- cutoffSD
  res$cutoffSD   <- out$cutoffSD
  out$highSDnai  <- which(out$SD > cutoffSD)
  res$highOD     <- out$highODnai  
  res$highSD     <- out$highSDnai
  NAimp          <- c(NAimp, out); rm(out)  
  
  if (bigOutput) {
    # Fully imputed data Xfi : scores and distances
    scoresfi <- (Xfi - matrix(rep(center, times = n), nrow = n, 
                              byrow = TRUE)) %*% loadings
    Fullimp  <- list(scoresfi = scoresfi)
    out      <- pca.distancesNew(res, Xfi, scoresfi, 
                                 XciSVD$rank, distprob)
    names(out)[1] <- "ODfi"
    names(out)[3] <- "SDfi"
    out$cutoffOD  <- cutoffOD
    out$highODfi  <- which(out$ODfi > cutoffOD)
    out$highSDfi  <- which(out$SDfi > cutoffSD)
    Fullimp  <- c(Fullimp, out); rm(out)
  }
  
  #########################################
  # Step 8: Unstandardization and residuals  
  #########################################
  
  # Compute residuals
  # C is for Centering:
  XOC <- sweep(XO, 2, center) # has NA's
  if (bigOutput) {
    XnaiC <- sweep(Xnai, 2, center)
    XciC  <- sweep(Xci, 2, center)
    XfiC  <- sweep(Xfi, 2, center)
  }
  
  # Compute standardized residuals of XO (with NA's): 
  if(k < XciSVD$rank){
    stdResid <- XOC - (scoresnai %*% t(loadings))
    res$residScale <- cellWise::estLocScale(stdResid, 
                                            nLocScale = DDCpars$nLocScale,
                                            type = "1stepM",
                                            precScale = DDCpars$precScale,
                                            center = FALSE)$scale
    res$stdResid <- sweep(stdResid, 2, res$residScale, "/")
    res$indcells <- which(abs(res$stdResid) >
                            sqrt(qchisq(DDCpars$tolProb, 1)))
    
    if (bigOutput) {
      # Compute standardized residuals of NA-imputed data:
      stdResidnai <- XnaiC - (scoresnai %*% t(loadings))
      NAimp$residScalenai <- cellWise::estLocScale(
        stdResidnai + 0*res$stdResid, # puts in NA's
        nLocScale = DDCpars$nLocScale, type = "1stepM",
        precScale = DDCpars$precScale, center = FALSE)$scale
      NAimp$stdResidnai <- sweep(stdResidnai, 2, NAimp$residScalenai, "/")
      NAimp$indcellsnai <- which(abs(NAimp$stdResidnai) >
                                   sqrt(qchisq(DDCpars$tolProb, 1)))
      
      # Compute standardized residuals of cell-imputed data:
      stdResidci <- XciC - (scoresci %*% t(loadings))
      Cellimp$residScaleci <- cellWise::estLocScale(
        stdResidci + 0*res$stdResid, # puts in NA's, 
        nLocScale = DDCpars$nLocScale, type = "1stepM",
        precScale = DDCpars$precScale, center = FALSE)$scale
      Cellimp$stdResidci <- sweep(stdResidci, 2, Cellimp$residScaleci, "/")
      Cellimp$indcellsci <- which(abs(Cellimp$stdResidci) >
                                    sqrt(qchisq(DDCpars$tolProb, 1)))
      
      # Compute standardized residuals of fully imputed data:
      stdResidfi <- XfiC - (scoresfi %*% t(loadings))
      Fullimp$residScalefi <- cellWise::estLocScale(
        stdResidfi + 0*res$stdResid, # puts in NA's
        nLocScale = DDCpars$nLocScale, type = "1stepM",
        precScale = DDCpars$precScale, center = FALSE)$scale
      Fullimp$stdResidfi <- sweep(stdResidfi, 2, Fullimp$residScalefi, "/")
      Fullimp$indcellsfi <- which(abs(Fullimp$stdResidfi) >
                                    sqrt(qchisq(DDCpars$tolProb, 1))) 
    }
  } else { # for k >= rank all residuals are zero, so we
    # cannot standardize them.
    res$stdResid <- 0*XOC # keeps the NA's
    res$residScale <- rep(0,d)
    res$indcells <- integer(0)
    if (bigOutput) {
      # residuals of NA-imputed data:
      NAimp$stdResidnai <- 0*XnaiC
      NAimp$residScalenai <- rep(0,d) 
      NAimp$indcellsnai <- integer(0)
      # residuals of cell-imputed data:
      Cellimp$stdResidci <- 0*XciC
      Cellimp$residScaleci <- rep(0,d)
      Cellimp$indcellsci <- integer(0)
      # residuals of fully imputed data:
      Fullimp$stdResidfi <- 0*XfiC
      Fullimp$residScalefi <- rep(0,d)
      Fullimp$indcellsfi <- integer(0)
    }
  }  
  
  ## put scales back:
  res$X.NAimp <- sweep(Xnai, 2, scaleX, "*")
  if (bigOutput) {
    Cellimp$Xci <- sweep(Xci, 2, scaleX, "*")  
    Fullimp$Xfi <- sweep(Xfi, 2, scaleX, "*")
  }
  
  ## update center:
  if(fixedCenter == TRUE){
    res$center  <- givenCenter
    res$X.NAimp <- sweep(res$X.NAimp, 2, givenCenter, "+")
    if (bigOutput) {
      Cellimp$Xci <- sweep(Cellimp$Xci, 2, givenCenter, "+")
      Fullimp$Xfi <- sweep(Fullimp$Xfi, 2, givenCenter, "+")      
    }
  } else {
    res$center <- center * scaleX
  }
  
  names(res$center)     <- colnames(checkOut$remX)
  names(res$scaleX)     <- colnames(checkOut$remX)
  names(res$residScale) <- colnames(checkOut$remX)
  if (bigOutput) {
    names(NAimp$residScalenai)  <- colnames(checkOut$remX)
    names(Cellimp$residScaleci) <- colnames(checkOut$remX)
    names(Fullimp$residScalefi) <- colnames(checkOut$remX)
  }
  
  # add remainder of output:
  if (bigOutput) {
    res$NAimp   <- NAimp
    res$Cellimp <- Cellimp   
    res$Fullimp <- Fullimp
  }
  return(res)
}


## AUXILIARY FUNCTIONS:

h.alpha.n <- function(alpha, n, p) {
  n2 <- (n + p + 1) %/% 2
  floor(2 * n2 - n + 2 * (n - n2) * alpha)
}

vecnorm <- function(x, p = 2) {
  sum(x ^ p) ^ (1 / p) 
}



genImpTable <- function(indimp, n, d){
  # gives indices of Imputed/Unimputed values per row,
  # based on the bivariate positions in indimp
  arrayimp <- arrayInd(indimp, c(n, d))
  imptab   <- vector("list", n)
  improws  <- seq_len(n)
  impcols  <- seq_len(d)
  for (i in improws) {   
    I  <- matrix(arrayimp[arrayimp[, 1] == i, ],ncol = 2)[, 2] 
    # Imputed variables
    U  <- setdiff(impcols, I) # Unmputed variables
    nI <- length(I)          # number of Imputed variables
    nU <- length(U)          # number of Unimputed variables 
    imptab[[i]] <- list(U = U, I = I, nU = nU, nI = nI)
  }
  return(imptab)
}


pca.distancesNew <- function(obj, data,myscores, r, crit = 0.99,
                             cutOD = NULL) { 
  n <- nrow(data)
  q <- ncol(myscores)
  smat <- diag(obj$eigenvalues,ncol = q)
  nk <- min(q, robustbase::rankMM(smat))
  if (nk < q) warning(paste("The smallest eigenvalue is ", 
                            obj$eigenvalues[q], 
                            " so the diagonal matrix of the eigenvalues",
                            " cannot be inverted!", sep = ""))
  SD <- sqrt(mahalanobis(as.matrix(myscores[, seq_len(nk)]), 
                         rep(0, nk), diag(obj$eigenvalues[seq_len(nk)],
                                          ncol = nk)))
  cutoffSD <- sqrt(qchisq(crit, obj$k))
  OD <- apply(data - matrix(rep(obj$center, times = n), nrow = n, 
                            byrow = TRUE) - 
                myscores %*% t(obj$loadings), 1, vecnorm)
  if (is.list(dimnames(myscores))) {
    names(OD) <- dimnames(myscores)[[1]]
  }
  if (is.null(cutOD)) {
    cutoffOD <- 0 # initializes
    if (obj$k != r) { # here we call critOD() :
      cutoffOD <- critOD(OD, crit = crit, umcd = TRUE, 
                         alpha = obj$alpha, classic = FALSE) }
  } else {cutoffOD <- cutOD }
  out <- list(OD = OD, cutoffOD = cutoffOD, SD = SD, cutoffSD = cutoffSD)
  return(out)
}


critOD <- function(OD, crit = 0.99, umcd = FALSE, alpha, classic = FALSE) 
{ # from rrcov:::crit.od
  OD <- OD ^ (2 / 3)
  if (classic) {
    t <- mean(OD)
    s <- sd(OD)
  }
  else if (umcd) {
    ms <- unimcd(OD, alpha)
    t <- ms$loc
    s <- ms$scale
  }
  else {
    t <- median(OD)
    s <- mad(OD)
  }
  cv <- (t + s * qnorm(crit)) ^ (3 / 2)
  cv
}


makeDirections <- function(data, ndirect, all = TRUE, 
                           fixedCenter = FALSE)
{
  uniran <- function(seed = 0) {
    seed <- floor(seed * 5761) + 999
    quot <- floor(seed/65536)
    seed <- floor(seed) - floor(quot * 65536)
    random <- seed/65536
    list(seed = seed, random = random)
  }
  randomset <- function(n, k, seed) {
    ranset <- vector(mode = "numeric", length = k)
    for (j in seq_len(k)) {
      r <- uniran(seed)
      seed <- r$seed
      num <- floor(r$random * n) + 1
      if (j > 1) {
        while (any(ranset == num)) {
          r <- uniran(seed)
          seed <- r$seed
          num <- floor(r$random * n) + 1
        }
      }
      ranset[j] <- num
    }
    ans <- list()
    ans$seed <- seed
    ans$ranset <- ranset
    ans
  }
  if (fixedCenter) {
    if (nrow(data) >= ndirect) { 
      # when n >= maxdir we draw a random sample from the rows:
      set.seed(0)
      sseed <- sample(1:nrow(data), size = ndirect,
                      replace = FALSE)
      B2 <- data[sseed, ]
    }
    else { # when n < maxdir we have to add more directions:
      ## In the main code we already took out the directions
      ## where vecnorm is zero first.
      data <- sweep(data, 1, apply(data, 1, function(x) sqrt(sum(x^2))), "/")
      # B2 <- data
      B2_add <- matrix(0, ndirect - nrow(data), ncol(data))
      seed <- 0
      r <- 1
      while (r <= ndirect - nrow(data)) {
        sseed <- randomset(nrow(data), 2, seed)
        seed <- sseed$seed
        B2_add[r, ] <- (data[sseed$ranset[1], ] + data[sseed$ranset[2], ])/2
        r <- r + 1
      }
      B2 <- rbind(data, B2_add)
    }
  } # ends fixedCenter == T   
  else { # no fixed center
    if (all) {
      cc <- utils::combn(nrow(data), 2)
      B2 <- data[cc[1, ], ] - data[cc[2, ], ]
    }
    else {
      n <- nrow(data)
      p <- ncol(data)
      r <- 1
      B2 <- matrix(0, ndirect, p)
      seed <- 0
      while (r <= ndirect) {
        sseed <- randomset(n, 2, seed)
        seed <- sseed$seed
        B2[r, ] <- data[sseed$ranset[1], ] - data[sseed$ranset[2],
        ]
        r <- r + 1
      }
    }
  }
  return(B2)
}


maxAngle <- function(Ptrue, Phat) {
  IPPI <- t(Ptrue) %*% Phat %*% t(Phat) %*% Ptrue
  lambdas_k <- eigen(IPPI, symmetric = TRUE)
  lambda_k_min <- min(lambdas_k$values)
  if (lambda_k_min > 1) lambda_k_min = 1
  angle <- acos((sqrt(lambda_k_min))) / (pi / 2) 
  return(angle)
}




truncPC = function (X, ncomp = NULL, scale = FALSE, center = TRUE, 
                    signflip = TRUE, via.svd = NULL, scores = FALSE) { 
  if (!is.numeric(X)) 
    stop(" X must be numeric.")
  Y <- as.matrix(X)
  if (!(length(dim(Y)) == 2)) 
    stop(" X is not a matrix.")
  n <- nrow(Y)
  d <- ncol(Y)
  if (n < 2) 
    stop(" The number of rows must be at least 2.")
  Y <- scale(Y, center = center, scale = scale)
  if (isTRUE(scale)) 
    scale <- attr(Y, "scaled:scale")
  if (isTRUE(center)) { center <- attr(Y, "scaled:center")
  } else { center = rep(0, d) }
  if (is.null(ncomp)) {
    SvdY <- try(svd::propack.svd(Y, neig = min(n, d)), silent = TRUE)
    if (inherits(SvdY, "try-error")) {
      SvdY <- svd(Y, nu = min(n, d), nv = min(n, d))
    }
  }
  else {
    if (!(ncomp >= 1)) 
      stop(" ncomp must be at least 1")
    ncomp <- min(ncomp, n, d)
    SvdY <- try(svd::propack.svd(Y, neig = ncomp), silent = TRUE)
    if (inherits(SvdY, "try-error")) {
      SvdY <- svd(Y, nu = ncomp, nv = ncomp)
    }
  }
  eigvals <- SvdY$d
  rank <- sum(eigvals > 1e-10)
  if (rank == 0) 
    stop(" The data has rank zero")
  eigvals <- (eigvals[seq_len(rank)])^2/(n - 1)
  loadings <- SvdY$v
  dim(loadings)
  loadings <- loadings[, seq_len(rank), drop = FALSE]
  if (signflip) {
    flipcolumn <- function(x) {
      if (x[which.max(abs(x))] < 0) {
        -x
      }
      else {
        x
      }
    }
    loadings <- apply(loadings, 2L, FUN = flipcolumn)
  }
  list(rank = rank, eigenvalues = eigvals, loadings = loadings, 
       scores = if (scores) Y %*% loadings, center = center, 
       scale = scale)
}


unimcd <- function(y, alpha = NULL, center = TRUE) {
  
  if (is.null(alpha)) {
    alpha <- 0.5
  }
  center <- as.numeric(center)
  res <- tryCatch(.Call('_cellWise_unimcd_cpp', y, alpha, center, PACKAGE = 'cellWise'),
                  "std::range_error" = function(e){conditionMessage(e)})
  return(list(loc = res$loc, scale = res$scale, weights = drop(res$weights)))
}
