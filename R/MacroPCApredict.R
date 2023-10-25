MacroPCApredict <- function(Xnew, InitialMacroPCA, MacroPCApars = NULL) { 
  if (is.null(MacroPCApars)) {
    MacroPCApars <- InitialMacroPCA$MacroPCApars
  }
  else {
    if (!is.list(MacroPCApars)) {
      stop("MacroPCApars must be a list")
    }
    InitialMacroPCA$MacroPCApars[names(MacroPCApars)] <- MacroPCApars
    MacroPCApars <- InitialMacroPCA$MacroPCApars
  }
  if (!"distprob" %in% names(MacroPCApars)) {
    MacroPCApars$distprob <- 0.99
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
  Xnew <- as.matrix(Xnew)
  n <- nrow(Xnew)
  d <- ncol(Xnew)
  dold <- ncol(InitialMacroPCA$remX)
  if(d != dold) stop(paste0(
    "Xnew should have ",dold," columns, that correspond\n",
    "  to those of InitialMacroPCA$remX"))
  oldcolnames <- colnames(InitialMacroPCA$remX)
  newcolnames <- colnames(Xnew)
  # print(oldcolnames)
  # print(newcolnames)
  if(!is.null(oldcolnames) & !is.null(newcolnames)){
    if(all.equal(newcolnames, oldcolnames) != TRUE) stop(
      paste0("Xnew should have the same column names",
             " as InitialMacroPCA$remX"))
  }
  InitialDDC  <- InitialMacroPCA$DDC
  DDCpars     <- MacroPCApars$DDCpars
  scaleX      <- InitialMacroPCA$scaleX
  distprob    <- MacroPCApars$distprob
  maxiter     <- MacroPCApars$maxiter
  tol         <- MacroPCApars$tol
  bigOutput   <- MacroPCApars$bigOutput
  resultDDC   <- DDCpredict(Xnew, InitialDDC, DDCpars)
  DDCpars     <- resultDDC$DDCpars
  DDCimp      <- sweep(resultDDC$Ximp, 2, scaleX, "/")
  XO          <- sweep(Xnew, 2, scaleX, "/"); rm(Xnew)
  Xfi         <- XO
  indcells    <- resultDDC$indcells
  indNA       <- which(is.na(XO))
  indimp      <- unique(c(indcells, indNA))
  Xfi[indimp] <- DDCimp[indimp] # impute flagged cells and NA's
  Xind        <- Xfi
  Xind[indimp] <- NA 
  Xind        <- is.na(Xind) # matrix with T in those cells
  k           <- InitialMacroPCA$k
  loadings    <- InitialMacroPCA$loadings
  center      <- InitialMacroPCA$center/scaleX
  eigenvalues <- InitialMacroPCA$eigenvalues
  if (any(Xind) & maxiter > 0) { # iterate, but keep center.
    diff <- 2 * tol
    It <- 0
    while (It < maxiter & diff > tol) {
      It <- It + 1
      Xfimis <- Xfi[Xind]
      XfiC <- sweep(Xfi, 2, center)           # centering
      Tr <- XfiC %*% loadings                 # scores
      Xfihat <- Tr %*% t(loadings)            # centered predictions
      Xfihat <- sweep(Xfihat, 2, center, "+") # uncentering
      Xfi[Xind] <- Xfihat[Xind]
      diff <- mean((Xfi[Xind] - Xfimis)^2)
    }
  }
  else {
    diff <- 0
    It <- 0
  }
  rownames(loadings) <- colnames(XO)
  colnames(loadings) <- paste0("PC", seq_len(ncol(loadings)))
  res <- list(MacroPCApars = MacroPCApars, DDC = resultDDC,
              scaleX = scaleX, k = k, loadings = loadings, 
              eigenvalues = eigenvalues, center = center, 
              It = It, diff = diff)
  Xnai <- XO
  Xnai[indNA]    <- Xfi[indNA]
  res$Xnew.NAimp <- sweep(Xnai, 2, scaleX, "*")
  XnaiC          <- sweep(Xnai, 2, center)
  scoresnai      <- XnaiC %*% loadings
  res$scores     <- scoresnai
  NAimp          <- list(scoresnai = scoresnai)
  cutoffOD       <- InitialMacroPCA$cutoffOD
  out <- pca.distancesNew(res, Xnai, scoresnai, ncol(Xnai), 
                          distprob, cutOD = cutoffOD)
  res$OD         <- out$OD
  out$cutoffOD   <- cutoffOD
  res$cutoffOD   <- out$cutoffOD
  res$SD         <- out$SD
  out$cutoffSD   <- InitialMacroPCA$cutoffSD
  res$cutoffSD   <- out$cutoffSD
  res$highOD     <- which(out$OD > out$cutoffOD)
  res$highSD     <- which(out$SD > out$cutoffSD)
  NAimp          <- c(NAimp, out)
  rm(out)
  XOC            <- sweep(XO, 2, center)
  stdResid       <- XOC - (scoresnai %*% t(loadings))
  res$residScale <- InitialMacroPCA$residScale
  res$stdResid   <- sweep(stdResid, 2, res$residScale, "/")
  cellcutoff     <- sqrt(qchisq(DDCpars$tolProb, 1))
  res$indcells   <- which(abs(res$stdResid) > cellcutoff)
  # res$Xnew.NAimp    <- sweep(Xnai, 2, scaleX, "*")
  res$center     <- center * scaleX # replaces earlier value
  names(res$scaleX)     <- colnames(Xnai)
  names(res$center)     <- colnames(Xnai)
  names(res$residScale) <- colnames(Xnai)
  #
  if (bigOutput) {
    stdResidnai <- XnaiC - (scoresnai %*% t(loadings))
    NAimp$residScalenai <- InitialMacroPCA$residScale
    NAimp$stdResidnai <- sweep(stdResidnai, 2, NAimp$residScalenai, 
                               "/")
    NAimp$indcellsnai <- which(abs(NAimp$stdResidnai) > cellcutoff)
    names(NAimp$residScalenai) <- colnames(Xnai)
    res$NAimp <- NAimp
    #
    Xci <- Xfi
    Xci[res$highOD] <- Xnai[res$highOD]
    XciC <- sweep(Xci, 2, center)
    scoresci <- XciC %*% loadings
    Cellimp <- list(scoresci = scoresci)
    out <- pca.distancesNew(res, Xci, scoresci, ncol(Xci), 
                            distprob, cutOD = cutoffOD)
    names(out)[1] <- "ODci"
    names(out)[3] <- "SDci"
    out$cutoffOD <- cutoffOD
    out$cutoffSD <- InitialMacroPCA$cutoffSD
    out$highODci <- which(out$ODci > out$cutoffOD)
    out$highSDci <- which(out$SDci > out$cutoffSD)
    Cellimp <- c(Cellimp, out)
    rm(out)
    stdResidci <- XciC - (scoresci %*% t(loadings))
    Cellimp$residScaleci <- InitialMacroPCA$Cellimp$residScaleci
    Cellimp$stdResidci <- sweep(stdResidci, 2, Cellimp$residScaleci, 
                                "/")
    Cellimp$indcellsci <- which(abs(Cellimp$stdResidci) > cellcutoff)
    Cellimp$Xci <- sweep(Xci, 2, scaleX, "*")
    names(Cellimp$residScaleci) <- colnames(Xnai)
    res$Cellimp <- Cellimp
    #
    XfiC          <- sweep(Xfi, 2, center)
    scoresfi      <- XfiC %*% loadings
    Fullimp       <- list(scoresfi = scoresfi)
    out <- pca.distancesNew(res, Xfi, scoresfi, ncol(Xfi), 
                            distprob, cutOD = cutoffOD)
    names(out)[1] <- "ODfi"
    names(out)[3] <- "SDfi"
    out$cutoffOD  <- cutoffOD
    out$cutoffSD  <- InitialMacroPCA$cutoffSD
    out$highODfi  <- which(out$ODfi > out$cutoffOD)
    out$highSDfi  <- which(out$SDfi > out$cutoffSD)
    Fullimp       <- c(Fullimp, out)
    rm(out)
    stdResidfi <- XfiC - (scoresfi %*% t(loadings))
    Fullimp$residScalefi <- InitialMacroPCA$Fullimp$residScalefi
    Fullimp$stdResidfi <- sweep(stdResidfi, 2, Fullimp$residScalefi, 
                                "/")
    Fullimp$indcellsfi <- which(abs(Fullimp$stdResidfi) > cellcutoff)
    Fullimp$Xfi <- sweep(Xfi, 2, scaleX, "*")
    names(Fullimp$residScalefi) <- colnames(Xnai)
    res$Fullimp <- Fullimp
  }
  return(res)
}