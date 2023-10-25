transfo = function(X, type = "YJ", robust = TRUE,
                   standardize = TRUE, 
                   quant = 0.99, nbsteps = 2, 
                   checkPars = list()) {
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  X <- as.matrix(X)
  dorig = ncol(X) # added
  colnamX = colnames(X)
  if(is.null(colnamX)) { colnamX = seq_len(dorig) }
  
  if (nbsteps < 1) {
    stop("nbsteps should be at least 1.")
  }
  if (!"coreOnly" %in% names(checkPars)) {
    checkPars$coreOnly <- FALSE
  }
  if (!"silent" %in% names(checkPars)) {
    checkPars$silent <- FALSE
  }
  if (!"numDiscrete" %in% names(checkPars)) {
    checkPars$numDiscrete <- 5
  }
  if (!"precScale" %in% names(checkPars)) {
    checkPars$precScale <- 1e-12
  }
  out <- NULL
  if (!checkPars$coreOnly) {
    out <- checkDataSet(X, fracNA = 1, numDiscrete = checkPars$numDiscrete, 
                        precScale = checkPars$precScale, 
                        silent = checkPars$silent)
    X <- out$remX
  }
  n <- nrow(X)
  d <- ncol(X) # from here on, only that many columns and entries
  Xt <- X
  lambdahats <- rep(1, d)
  ttypes <- rep(NA, d)
  objective <- rep(NA, d)
  weights <- matrix(0, n, d)
  muhat <- locx <- rep(0, d)
  sigmahat <- scalex <- rep(1, d)
  lambdarange <- c(-4, 6)
  for (j in seq_len(d)) { # j=1
    x <- X[, j]
    goodInds <- which(!is.na(x))
    x <- x[goodInds]
    ttype <- getTtype(x, type, j)
    ttypes[j] <- ttype
    if (robust == TRUE) {
      order.x <- order(x)
      xsort <- x[order.x]
      if (ttype == "BC" || ttype == "bestObj") {
        lambdarangetemp <- lambdarange
        converged <- FALSE
        while (!converged) {
          est.out.BC <- RewML_BC(xsort, 
                                 lambdarange = lambdarangetemp, 
                                 standardize = standardize, 
                                 init = "BCr", 
                                 quant = quant, nbsteps = nbsteps)
          converged <- min(abs(est.out.BC$lambdahat.rew - 
                                 lambdarangetemp)) > diff(lambdarangetemp) * 
            0.05
          lambdarangetemp <- 1 + (lambdarangetemp - 1) * 
            (1 + (abs(est.out.BC$lambdahat.rew - lambdarangetemp) == 
                    min(abs(est.out.BC$lambdahat.rew - lambdarangetemp))) + 
               0)
        }
        lambdahats[j] <- est.out.BC$lambdahat.rew
        objective[j] <- est.out.BC$critval.rew
        Xt[goodInds[order.x], j] <- est.out.BC$yt.rew
        weights[goodInds[order.x], j] <- est.out.BC$weights
        muhat[j] <- est.out.BC$muhat
        sigmahat[j] <- est.out.BC$sigmahat
        locx[j] <- est.out.BC$locx
        scalex[j] <- est.out.BC$scalex
      }
      if (ttype == "YJ" || ttype == "bestObj") {
        lambdarangetemp <- lambdarange
        converged <- FALSE
        while (!converged) {
          est.out.YJ <- RewML_YJ(xsort,
                                 lambdarange = lambdarangetemp, 
                                 standardize = standardize, 
                                 init = "YJr", 
                                 quant = quant, nbsteps = nbsteps)
          converged <- min(abs(est.out.YJ$lambdahat.rew - 
                                 lambdarangetemp)) > diff(lambdarangetemp) * 
            0.05
          lambdarangetemp <- 1 + (lambdarangetemp - 1) * 
            (1 + (abs(est.out.YJ$lambdahat.rew - lambdarangetemp) == 
                    min(abs(est.out.YJ$lambdahat.rew - lambdarangetemp))) + 
               0)
        }
        lambdahats[j] <- est.out.YJ$lambdahat.rew
        objective[j] <- est.out.YJ$critval.rew
        Xt[goodInds[order.x], j] <- est.out.YJ$yt.rew
        weights[goodInds[order.x], j] <- est.out.YJ$weights
        muhat[j] <- est.out.YJ$muhat
        sigmahat[j] <- est.out.YJ$sigmahat
        locx[j] <- est.out.YJ$locx
        scalex[j] <- est.out.YJ$scalex
      }
      if (ttype == "bestObj") {
        if (est.out.BC$critval.rew < est.out.YJ$critval.rew) {
          lambdahats[j] <- est.out.BC$lambdahat.rew
          objective[j] <- est.out.BC$critval.rew
          Xt[goodInds[order.x], j] <- est.out.BC$yt.rew
          weights[goodInds[order.x], j] <- est.out.BC$weights
          ttypes[j] <- "BC"
          muhat[j] <- est.out.BC$muhat
          sigmahat[j] <- est.out.BC$sigmahat
          locx[j] <- est.out.BC$locx
          scalex[j] <- est.out.BC$scalex
          ttypes[j] <- "BC"
        }
        else {
          ttypes[j] <- "YJ"
        }
      }
    }
    else { # robust = F
      if (ttype == "BC" || ttype == "bestObj") {
        lambdarangetemp <- lambdarange
        converged <- FALSE
        while (!converged) {
          est.out.BC <- estML(x, type = "BC", 
                              lambdarange = lambdarangetemp, 
                              standardize = standardize)
          converged <- min(abs(est.out.BC$lambda - lambdarangetemp)) > 
            diff(lambdarangetemp) * 0.05
          lambdarangetemp <- 1 + (lambdarangetemp - 1) * 
            (1 + (abs(est.out.BC$lambda - lambdarangetemp) == 
                    min(abs(est.out.BC$lambda - lambdarangetemp))) + 
               0)
        }
        lambdahats[j] <- est.out.BC$lambda
        objective[j] <- est.out.BC$objective
        locx[j] <- est.out.BC$locx
        scalex[j] <- est.out.BC$scalex
        Xt[goodInds, j] <- est.out.BC$xt
        weights[goodInds, j] <- est.out.BC$weights
      }
      if (ttype == "YJ" || ttype == "bestObj") {
        lambdarangetemp <- lambdarange
        converged <- FALSE
        while (!converged) {
          est.out.YJ <- estML(x, type = "YJ", 
                              lambdarange = lambdarangetemp, 
                              standardize = standardize)
          converged <- min(abs(est.out.YJ$lambda - lambdarangetemp)) > 
            diff(lambdarangetemp) * 0.05
          lambdarangetemp <- 1 + (lambdarangetemp - 1) * 
            (1 + (abs(est.out.YJ$lambda - lambdarangetemp) == 
                    min(abs(est.out.YJ$lambda - lambdarangetemp))) + 
               0)
        }
        lambdahats[j] <- est.out.YJ$lambda
        objective[j] <- est.out.YJ$objective
        locx[j] <- est.out.YJ$locx
        scalex[j] <- est.out.YJ$scalex
        Xt[goodInds, j] <- est.out.YJ$xt
        weights[goodInds, j] <- est.out.YJ$weights
      }
      if (ttype == "bestObj") {
        if (est.out.BC$objective > est.out.YJ$objective) {
          lambdahats[j] <- est.out.BC$lambda
          objective[j] <- est.out.BC$objective
          locx[j] <- est.out.BC$locx
          scalex[j] <- est.out.BC$scalex
          Xt[goodInds, j] <- est.out.BC$xt
          weights[goodInds, j] <- est.out.BC$weights
          ttypes[j] <- "BC"
        }
        else {
          ttypes[j] <- "YJ"
        }
      }
      muhat[j] <- mean(Xt[, j], na.rm = TRUE)
      sigmahat[j] <- sd(Xt[, j], na.rm = TRUE)
    }
  }
  # Zt <- scale(Xt, center = muhat, scale = sigmahat)
  Y <- if (standardize) { scale(Xt, center = muhat, scale = sigmahat)
  }
  else {
    Xt
  }
  # must add dorig to output list
  return(c(list(lambdahats = lambdahats, objective = objective, 
                Y = Y, weights = weights, ttypes = ttypes, muhat = muhat, 
                sigmahat = sigmahat, locx = locx, scalex = scalex, 
                standardize = standardize), out, list(dorig = dorig,
                                                      colnamX = colnamX) ) )
}


transfo_newdata = function(Xnew, transfo.out) {
  if (is.vector(Xnew)) {
    Xnew <- matrix(Xnew, ncol = 1)
  }
  Xnew <- as.matrix(Xnew)
  d = ncol(Xnew)
  dorig = transfo.out$dorig
  if(d != dorig) stop(paste0(
    "Xnew should have ",dorig," columns like the original data."))
  colnamnew = colnames(Xnew)
  if(is.null(colnamnew)) { colnamnew = seq_len(d) }
  if(all.equal(colnamnew, transfo.out$colnamX) != TRUE) stop(
    "Xnew should have the same column names as the original data.")
  colInA <- seq_len(d)
  if (!is.null(transfo.out$out$colInAnalysis)) {
    colInA <- transfo.out$out$colInAnalysis
  }
  standardize <- transfo.out$standardize
  locx <- transfo.out$locx
  scalex <- transfo.out$scalex
  muhat <- transfo.out$muhat
  sigmahat <- transfo.out$sigmahat
  ttypes <- transfo.out$ttypes
  lambdas <- transfo.out$lambdahats
  Ynew <- Xnew
  for (j in colInA) {
    xnewt <- Xnew[, colInA[j]] # actual column
    if (ttypes[j] == "YJ") {
      ynewt <- xnewt
      if (standardize) {
        ynewt <- scale(ynewt, locx[j], scalex[j])
      }
      ynewt <- YJ(ynewt, lambdas[j], stdToo = FALSE)$yt
      if (standardize) {
        ynewt <- scale(ynewt, muhat[j], sigmahat[j])
      }
      Ynew[, j] <- ynewt
    }
    else if (ttypes[j] == "BC") {
      ynewt <- xnewt
      if (standardize) {
        ynewt <- ynewt/scalex[j]
      }
      if(min(ynewt) <= 0) stop(paste0(
        "Column ",j," has some value(s) <= 0, but in the original data\n",
        "it was transformed by Box-Cox. You have to either make all\n",
        "values in this column strictly positive, or apply transfo()\n",
        "again to the original data but with type = \"YJ\"."))
      ynewt <- BC(ynewt, lambdas[j], stdToo = FALSE)$yt
      if (standardize) {
        ynewt <- scale(ynewt, muhat[j], sigmahat[j])
      }
      Ynew[, j] <- ynewt
    }
    else {
      stop(paste0("Invalid transformation type ",ttypes[j],
                  " for column ",j))
    }
  }
  return(Ynew)
}


transfo_transformback = function (Ynew, transfo.out) {
  if (is.vector(Ynew)) {
    Ynew <- matrix(Ynew, ncol = 1)
  }
  Ynew <- as.matrix(Ynew)
  d = ncol(Ynew)
  dorig = transfo.out$dorig
  if(d != transfo.out$dorig) stop(paste0(
    "Ynew should have ",dorig," columns like the original ",
    "transformed data."))
  colnamnew = colnames(Ynew)
  if(is.null(colnamnew)) { colnamnew = seq_len(d) }
  if(all.equal(colnamnew, transfo.out$colnamX) != TRUE) stop(paste0(
    "Ynew should have the same column names as the original ",
    "transformed data."))
  standardize <- transfo.out$standardize
  locx <- transfo.out$locx
  scalex <- transfo.out$scalex
  ttypes <- transfo.out$ttypes
  lambdas <- transfo.out$lambdahats
  muhat <- transfo.out$muhat
  sigmahat <- transfo.out$sigmahat
  colInAnalysis <- seq_len(ncol(Ynew))
  # if (!is.null(transfo.out$out$colInAnalysis)) {
  #   colInAnalysis <- transfo.out$out$colInAnalysis
  # }
  Xnew <- Ynew # initialization
  for (j in colInAnalysis) {
    # Shouldn't we step over j in seq_len(length(colInAnalysis)) ?
    ynewt <- Ynew[, j]
    if (ttypes[j] == "YJ") {
      xnewt <- ynewt
      stan = " "
      if (standardize) {
        xnewt <- xnewt * sigmahat[j] + muhat[j]
        stan = " destandardized "
      }
      if(lambdas[j] > 2){
        lowerb = -1/abs(lambdas[j] - 2)
        newlowerb = 0.95*lowerb
        if(min(xnewt) < lowerb) warning(paste0(
          "The lowest",stan,"value in column ",j," was ",min(xnewt),
          " .\n","This is below the expected lower bound ",lowerb,
          " .\n","Such low",stan,"values were put equal to ",
          newlowerb,"\nso they could be transformed back."))
        xnewt = pmax(xnewt, newlowerb)
      }
      if(lambdas[j] < 0){
        upperb = 1/abs(lambdas[j])
        newupperb = 0.95*upperb
        if(max(xnewt) > upperb) warning(paste0(
          "The highest",stan,"value in column ",j," was ",max(xnewt),
          " .\n","This is above the expected upper bound ",upperb,
          " .\n","Such high",stan,"values were put equal to ",
          newupperb,"\nso they could be transformed back."))
        xnewt = pmin(xnewt, newupperb)
      }
      xnewt <- iYJ(xnewt, lambdas[j], stdToo = FALSE)$yt
      if (standardize) {
        xnewt <- xnewt * scalex[j] + locx[j]
      }
      Xnew[, j] <- xnewt
    }
    else if (ttypes[j] == "BC") {
      xnewt <- ynewt
      stan = " "
      if (standardize) {
        xnewt <- xnewt * sigmahat[j] + muhat[j]
        stan = " destandardized "
      }
      if(lambdas[j] > 0){
        lowerb = -1/abs(lambdas[j])
        newlowerb = 0.95*lowerb
        if(min(xnewt) < lowerb) warning(paste0(
          "The lowest",stan,"value in column ",j," was ",min(xnewt),
          " .\n","This is below the expected lower bound ",lowerb,
          " .\n","Such low",stan,"values were put equal to ",
          newlowerb,"\nso they could be transformed back."))
        xnewt = pmax(xnewt, newlowerb)  
      }
      if(lambdas[j] < 0){
        upperb = 1/abs(lambdas[j])
        newupperb = 0.95*upperb
        if(max(xnewt) > upperb) warning(paste0(
          "The highest",stan,"value in column ",j," was ",max(xnewt),
          " .\n","This is above the expected upper bound ",upperb,
          " .\n","Such high",stan,"values were put equal to ",
          newupperb,"\nso they could be transformed back."))
        xnewt = pmin(xnewt, newupperb)
      }
      xnewt <- iBC(xnewt, lambdas[j], stdToo = FALSE)$yt
      if (standardize) {
        xnewt <- xnewt * scalex[j]
      }
      Xnew[, j] <- xnewt
    }
    else {
      stop(paste0("Invalid transformation type ",ttypes[j],
                  " for column ",j))
    }
  }
  return(Xnew)
}


getTtype <- function(x, type, j=NULL) {
  # Gets the transformation type (BC, YJ, or bestObj) and checks 
  # whether the domain of the variable is compatible with that type.
  # args:
  #   x    : vector with the values of the variable
  #   type : one of "BC", "YJ", "bestObj"
  #   j    : number of the variable
  # 
  if (!type %in% c("BC", "YJ", "bestObj")) {
    stop(" type should be 'BC' or 'YJ' or 'bestObj'.")
  }
  minx <- min(x)
  if (type == "BC") {
    if (minx <= 0) {
      stop(paste("The minimum of variable ", j,
                 " is not strictly positive.\n",
                 "Please choose a type other than BC.",sep = ""))
    } else {
      ttype <- "BC"
    }
  } else if (type == "YJ") {
    ttype <- "YJ"
  } else if (type == "bestObj") {
    if (minx <= 0) {
      ttype <- "YJ"
    } else {
      ttype <- "bestObj"
    }
  }
  return(ttype)
}


RewML_BC <- function(x,
                     lambdarange = NULL,
                     standardize = TRUE,
                     init = "BCr",
                     quant = 0.99,
                     nbsteps = 2) {
  # The reweighted ML estimator for the Box-Cox transformation parameter.
  #
  # args: 
  #   x:              vector of _sorted_ observations
  #   lambdarange:        grid of lambda values. If NULL, a grid between 
  #                   -2 and 4 is chosen
  #   init:           initial estimator. should be "BCr" or "BC"
  #   quant:          quantile for determining the weights in the 
  #                   reweighting step
  #   nbsteps:        number of reweighting steps
  #
  if (!init %in% c("BC", "BCr")) {
    stop("init should be either 'BC' or 'BCr'")
  }
  x <- na.omit(x)
  
  if (is.unsorted(x)) {
    x.order <- order(x, decreasing = FALSE)
    x       <- x[x.order]
  } else {
    x.order <- NULL
  }
  scalex <- 1
  if (standardize) {
    scalex <- median(x)
    x <- x / scalex # so median is 1
    # x <- exp(log(x)/mad(log(x)) ) # median stays 1
    #  if(min(x) < 1e-6) { x <- x + 1e-6 }
  }
  
  # Range of lambda over which to optimize:
  if (is.null(lambdarange)) { lambdarange <- c(-4, 6) }
  if (init == "BCr") {
    tempfunc <- function(lambdatemp) {
      getCritval(x, lambdatemp, tfunc = BCr)
    }
    opt.out <- optimize(tempfunc, range(lambdarange))
    lambdahat.raw <- opt.out$minimum
    BCr.out.raw   <- BCr(x, lambdahat.raw)
    yt.raw        <- BCr.out.raw$yt
    zt.raw        <- BCr.out.raw$zt
  } else { 
    tempfunc <- function(lambdatemp) {
      getCritval(x, lambdatemp, tfunc = BC)
    }
    opt.out <- optimize(tempfunc, range(lambdarange))
    lambdahat.raw <- opt.out$minimum
    BC.out.raw    <- BC(x, lambdahat.raw)
    yt.raw        <- BC.out.raw$yt
    zt.raw        <- BC.out.raw$zt
  }
  rew.out <- reweightBCr(x = x, zt.raw = zt.raw,
                         lambdahat.raw = lambdahat.raw, 
                         lambdarange = lambdarange, quant = quant,
                         nbsteps = nbsteps,
                         standardize = standardize)
  BC.out.rew    <- rew.out$BC.out.rew
  yt.rew        <- BC.out.rew$yt
  zt.rew        <- BC.out.rew$zt
  critval.rew   <- rew.out$critval.rew
  lambdahat.rew <- rew.out$lambdahat.rew
  weights       <- rew.out$wgts
  muhat         <- mean(yt.rew[weights == 1])
  sigmahat      <- sd(yt.rew[weights == 1])
  
  if (!is.null(x.order)) {
    yt.raw[x.order] <- yt.raw
    zt.raw[x.order] <- zt.raw
    yt.rew[x.order] <- yt.rew
    zt.rew[x.order] <- zt.rew
    weights[x.order] <- weights
  }
  return(list(lambdahat.raw = lambdahat.raw,
              yt.raw = yt.raw,
              zt.raw = zt.raw,
              weights = weights,
              lambdahat.rew = lambdahat.rew,
              yt.rew = yt.rew,
              zt.rew = zt.rew,
              critval.rew = critval.rew, 
              muhat = muhat,
              sigmahat = sigmahat,
              locx = 0,
              scalex = scalex))
}


RewML_YJ <- function(x, 
                     lambdarange = NULL, 
                     quant = 0.99, 
                     init = "YJr",
                     standardize = TRUE, 
                     nbsteps = 2) {
  # The reweighted ML estimator for the Yeo-Johnson transformation parameter.
  #
  # args: 
  #   x             : vector of _sorted_ observations
  #   lambdarange       : grid of lambda values. If NULL, a grid between 
  #                   -2 and 4 is chosen
  #   init          : initial estimator. should be "YJr" or "YJ"
  #   quant         : quantile for determining the weights in the
  #                   reweighting step
  #   nbsteps       : number of reweighting steps
  #
  if (!init %in% c("YJ", "YJr")) {
    stop("init should be either 'YJ' or 'YJr'")
  }
  
  x <- na.omit(x)
  
  if(is.unsorted(x)) {
    x.order <- order(x, decreasing = FALSE)
    x       <- x[x.order]
  } else {
    x.order <- NULL
  }
  
  locx <- 0
  scalex <- 1
  if (standardize) {
    locx <- median(x)
    scalex <- mad(x)
    x <- (x - locx) / scalex
  }
  
  # Range of lambda over which to optimize:
  if (is.null(lambdarange)) { lambdarange <- c(-4,6) }
  if (init == "YJr") {
    tempfunc <- function(lambdatemp) {
      getCritval(x, lambdatemp, tfunc = YJr)
    }
    opt.out <- optimize(tempfunc, range(lambdarange))
    lambdahat.raw <- opt.out$minimum
    YJr.out.raw   <- YJr(x, lambdahat.raw)
    yt.raw        <- YJr.out.raw$yt
    zt.raw        <- YJr.out.raw$zt
  } else { 
    tempfunc <- function(lambdatemp) {
      getCritval(x, lambdatemp, tfunc = YJ)
    }
    opt.out <- optimize(tempfunc, range(lambdarange))
    lambdahat.raw <- opt.out$minimum
    YJr.out.raw   <- YJ(x, lambdahat.raw)
    yt.raw        <- YJr.out.raw$yt
    zt.raw        <- YJr.out.raw$zt
  }
  # Reweighting:
  rew.out <- reweightYJr(x, zt.raw, lambdahat.raw, lambdarange, 
                         quant = quant, nbsteps = nbsteps)
  weights       <- rew.out$wgts
  YJ.out.rew    <- rew.out$YJ.out.rew
  yt.rew        <- YJ.out.rew$yt
  zt.rew        <- YJ.out.rew$zt
  critval.rew   <- rew.out$critval.rew
  lambdahat.rew <- rew.out$lambdahat.rew
  muhat         <- mean(yt.rew[weights == 1])
  sigmahat      <- sd(yt.rew[weights == 1])
  
  if (!is.null(x.order)) {
    yt.raw[x.order] <- yt.raw
    zt.raw[x.order] <- zt.raw
    yt.rew[x.order] <- yt.rew
    zt.rew[x.order] <- zt.rew
    weights[x.order] <- weights
  }
  
  return(list(lambdarange = lambdarange,
              lambdahat.raw = lambdahat.raw,
              yt.raw = yt.raw,
              zt.raw = zt.raw,
              weights = weights,
              critval.rew = critval.rew,
              lambdahat.rew = lambdahat.rew,
              yt.rew = yt.rew,
              zt.rew = zt.rew, 
              muhat = muhat,
              sigmahat = sigmahat,
              locx = locx,
              scalex = scalex))
}


reweightYJr <- function(x, 
                        zt.raw, 
                        lambdahat.raw,
                        lambdarange,
                        quant = 0.99, 
                        nbsteps = 2) {
  # Function for reweighted maximum likelihood, based on an initial estimate.
  # args: 
  #   x             : vector of sorted original observations
  #   zt.raw        : vector of sorted poststandardized transformed data 
  #                   (with initial lambdahat)
  #   lambdahat.raw : initial estimate for transformation parameter lambda
  #   lambdarange   : range of lambda values.
  #   quant         : quantile for determining the weights in the 
  #                   reweighting step
  #   nbsteps       : number of reweighting steps
  #
  if (is.null(lambdarange)) { lambdarange = c(-4,6) }
  # initial reweighting:
  wgts    <- abs(zt.raw) <= sqrt(qchisq(quant,1))
  rewinds <- which(wgts == 1)
  x.rew   <- x[rewinds]
  x.all   <- x
  # YJ transform with weighted ML:
  for (k in seq_len(nbsteps)) {
    lambdahat.rew <- estML(x = x.rew, lambdarange = lambdarange, 
                           type = "YJ", standardize = FALSE)$lambda
    YJ.out.rew    <- YJ(x.all, lambdahat.rew)
    wgts    <- abs(YJ.out.rew$zt) <= sqrt(qchisq(quant,1))
    rewinds <- which(wgts == 1)
    x.rew   <- x[rewinds]
  }
  critval.rew <- getCritval(x.all, lambdahat.rew, tfunc = YJ,
                            quant=quant)
  return(list(wgts = wgts,
              critval.rew  = critval.rew,
              lambdahat.rew = lambdahat.rew,
              YJ.out.rew = YJ.out.rew))
}


getChangepointYJr <- function(y, lambda, fac = 1.5, eps = 1e-5) {
  # Function to caculate the "changepoint" for the rectified
  # Yeo-Johnson transformation
  # 
  quarts = localFivenum(y)[2:4]
  if (lambda < 1){
    chg = YJ(quarts[3], lambda, stdToo = FALSE)$yt * 1.5
  } else if (lambda > 1){
    chg = YJ(quarts[1], lambda, stdToo = FALSE)$yt * 1.5
  }
  if (lambda < 0) {
    chg <- min(chg, abs(1 / lambda) - eps)
  } else if (lambda > 2) {
    chg <- max(chg, (1 / (2 - lambda)) + eps)
  }
  chg <- iYJ(chg, lambda, stdToo = FALSE)$yt
  chg <- min(max(chg, y[1]), y[length(y)])
  return(chg)
}


getChangepointBCr <- function(y, lambda, fac = 1.5, eps = 1e-5) {
  # Function to caculate the "changepoint" for the rectified
  # Box Cox transformation
  # 
  quarts = localFivenum(y)[2:4]
  if (lambda < 1){
    chg = BC(quarts[3], lambda, stdToo = FALSE)$yt * 1.5
  } else if (lambda > 1){
    chg = BC(quarts[1], lambda, stdToo = FALSE)$yt * 1.5
  }
  #correct for image of BC
  if (lambda < 0) {
    chg <- min(chg, abs(1 / lambda) - eps)
  } else  if (lambda > 0) {
    chg <- max(chg, -(1 / lambda) + eps)
  }
  chg <- iBC(chg, lambda, stdToo = FALSE)$yt
  chg <- min(max(chg, y[1]), y[length(y)])
  return(chg)
}


reweightBCr <- function(x, zt.raw, lambdahat.raw, lambdarange, quant = 0.99,
                        nbsteps = 2, standardize = TRUE) {
  # function for reweighted maximum likelihood, based on an initial estimate.
  # args: 
  #   x              : vector of sorted original observations
  #   zt.raw         : vector of sorted poststandardized transformed data 
  #                    (from the initial lambdahat)
  #   lambdahat.raw  : initial estimate for transformation parameter lambda
  #   lambdarange    : range of lambda values.
  #   quant          : quantile for determining the weights in the 
  #                    reweighting step
  #   nbsteps        : number of reweighting steps
  #
  if (is.null(lambdarange)) { lambdarange = c(-4,6) }
  # Initial weights:
  weights <- abs(zt.raw) <= sqrt(qchisq(quant,1))
  x.all   <- x
  x.rew   <- x[which(weights == 1)]
  # 
  # BC transform with weighted ML
  #
  for (k in seq_len(nbsteps)) {
    lambdahat.rew <- estML(x = x.rew, lambdarange = lambdarange, 
                           type = "BC", standardize = FALSE)$lambda
    BC.out.rew <- BC(x.all, lambdahat.rew)
    wgts       <- abs(BC.out.rew$zt) <= sqrt(qchisq(quant,1))
    rewinds    <- which(wgts == 1)
    x.rew      <- x[rewinds]
  }
  critval.rew <- getCritval(x.all, lambdahat.rew, tfunc = BC,quant = quant)
  return(list(BC.out.rew = BC.out.rew,
              lambdahat.rew = lambdahat.rew, 
              critval.rew = critval.rew,
              wgts = wgts))
}



estML <- function(x, lambdarange = NULL, type = "BC",
                  standardize = TRUE) {
  # Computes the ML estimator of the parameter lambda in the BC 
  # and YJ transformation, for a vector x (corresponding to a
  # column of the input data X).
  # Assumes x has no NA's, as they were already taken out in the main. 
  #
  n <- length(x)
  if (is.null(lambdarange)) { lambdarange <- c(-4, 6) }
  if (type == "YJ") { # Yeo-Johnson
    locx <- 0
    scalex <- 1
    if (standardize) {
      locx <- mean(x)
      scalex <- sd(x)
      x <- (x - locx) / scalex
    }
    tempfunc <- function(lambdatemp) {
      xyj    <- YJ(x, lambdatemp, stdToo = FALSE)$yt
      mu     <- mean(xyj)
      sigma2 <- mean((xyj - mu)^2)
      result <- -0.5 * n * log(2 * pi) -
        0.5 * n * log(sigma2) -
        0.5/sigma2 * sum( (xyj - mu)^2) +
        (lambdatemp - 1) * sum(sign(x) * log(1 + abs(x)))
      return(result)
    }
    opt.out <- optimize(tempfunc, range(lambdarange), maximum = TRUE)
    objective = tempfunc(opt.out$maximum) 
    xt <- YJ(x, opt.out$maximum, stdToo = FALSE)$yt
    zt <- scale(xt)
  } else {# Box-Cox
    locx <- 0
    scalex <- 1
    if (standardize) {
      scalex <- median(x)
      x <- x / scalex
      #x <- exp( log(x)/sd(log(x)) ) # median stays 1
      #if(min(x) < 1e-6) { x <- x + 1e-6 }
      
    }
    tempfunc <- function(lambdatemp) {
      xyj    <- BC(x, lambdatemp, stdToo = FALSE)$yt
      mu     <- mean(xyj)
      sigma2 <- mean((xyj - mu)^2)
      result <- -n/2 * log(sigma2) + (lambdatemp - 1) * sum(log(x))
      return(result)
    }
    opt.out <- optimize(tempfunc, range(lambdarange), maximum = TRUE)
    objective = tempfunc(opt.out$maximum) 
    xt <- BC(x, opt.out$maximum, stdToo = FALSE)$yt
    zt <- scale(xt)
  }
  return(list(lambda = opt.out$maximum,
              lambdarange = lambdarange, objective = objective,
              xt = xt, zt = zt, weights = rep(1, n),
              locx = locx, scalex = scalex))
}


getCritval <- function(x, lambda, tfunc = YJr, chg2=NULL,  quant=0.99) {
  # Calculates the value of the robust fitting criterion.
  # Arguments:
  # x      : the univariate data
  # lambda : the value of lambda
  # tfunc  : either BC, BCr, YJ, YJr
  # chg2   : 2 changepoints: one for lambda >= 1 and one if lambda < 1
  #          so chg[1] < chg[2]
  #          If NULL, uses Q1 for lambda > 1 and Q3 for lambda > 1
  #
  if (is.null(chg2)) { chg = NULL } else {
    chg <- if (lambda >= 1) {chg2[1]} else {chg2[2]}
  }
  robnormality(tfunc(x,lambda,chg,stdToo = FALSE)$yt, b = 0.5)
}


BC <- function(y, lambda, chg = NULL, stdToo = TRUE) {
  # The Box-Cox transformation.
  #
  if (lambda == 0) {
    yt <- log(y)
  } else {
    yt <- (y^lambda - 1) / lambda
  }
  if (stdToo) {
    if (length(y) > 1) {
      locScale <- estLocScale(matrix(yt, ncol = 1), 
                              type = "hubhub")
      zt <- (yt - locScale$loc) / locScale$scale
    } else { 
      zt <- yt
    }
  } else {
    zt <- NULL
  }
  return(list(yt = yt, zt =zt))
}


iBC <- function(y, lambda, chg = NULL, stdToo = TRUE) {
  # Inverse of the Box-Cox transformation.
  #
  if (lambda == 0) {
    yt <- exp(y)
  } else {
    yt <- (1 + y * lambda)^(1/lambda)
  }
  
  if (stdToo) {
    zt = (yt - median(yt)) / mad(yt)
  } else {
    zt = NULL
  }
  return(list(yt = yt, zt = zt))
}


YJ <- function(y, lambda, chg = NULL, stdToo = TRUE) {
  # The Yeo-Johnson transformation.
  #
  indlow  <- which(y < 0)
  indhigh <- which(y >= 0)
  if (lambda != 0) {
    y[indhigh] = ((1 + y[indhigh])^(lambda) - 1) / lambda
  } else { 
    y[indhigh] = log(1 + y[indhigh])
  }
  if (lambda != 2) {
    y[indlow] = -((1 - y[indlow])^(2 - lambda) - 1) / (2 - lambda)
  } else {
    y[indlow] = -log(1 - y[indlow])
  }
  if (stdToo) {
    if (length(y) > 1) {
      locScale <- estLocScale(matrix(y, ncol = 1), 
                              type = "hubhub")
      zt <- (y - locScale$loc) / locScale$scale
    } else {
      zt <- y
    }
  } else {
    zt <- NULL
  }
  return(list(yt = y, zt = zt))
}


iYJ <- function(y, lambda, stdToo = TRUE) {
  # Inverse of the Yeo-Johnson transformation.
  #
  indlow  <- which(y < 0)
  indhigh <- which(y >= 0)
  if (lambda != 0) {
    y[indhigh] = (1 + lambda * y[indhigh])^(1 / lambda) - 1
  } else { 
    y[indhigh] = exp(y[indhigh]) - 1
  }
  if (lambda != 2) {
    y[indlow] = -((1 + (lambda - 2) * y[indlow])^(1/(2-lambda))) + 1
  } else {
    y[indlow] = 1 - exp(-y[indlow])
  }
  if (stdToo) {
    zt = (y - median(y)) / mad(y)
  } else {
    zt = NULL
  }
  return(list(yt = y, zt = zt))
}


YJprimex <- function(x, lambda) {
  # partial derivative of YJ with respect to the argument x.
  (1 + abs(x))^(sign(x) * (lambda - 1))
}


estMTL <- function(x, alpha = 0.95,
                   lambdagrid = NULL, type = "YJ",
                   standardize = TRUE) {
  # The Maximum Trimmed Likelihood estimator.
  # 
  if (is.null(lambdagrid)) {
    lambdagrid <- ((0:200) - 100)/50
  }
  result  <- rep(0, length(lambdagrid))
  n       <- length(x)
  h       <- ceiling(alpha * n)
  if (standardize) {
    if (type == "YJ") {
      x <- (x - median(x))/mad(x)
    }
    if (type == "BC") {
      x <- x / median(x)
    }
  }
  x_order <- order(x) # ranks of the data
  if (type == "YJ") { # Yeo-Johnson
    for (i in seq_len(length(lambdagrid))) { # loop over lambdagrid
      lambda <- lambdagrid[i]
      xyj <- YJ(x, lambda)$yt
      nbSubsets <- n - h + 1
      sds   <- rep(1, nbSubsets)
      means <- rep(0, nbSubsets)
      for (j in seq_len(nbSubsets)) { # loop over subsets
        datatemp <- xyj[x_order[j:(j + h - 1)]]
        # subset of ordered xyj
        means[j] <- mean(datatemp)
        sds[j]   <- sd(datatemp)
      } # ends loop over subsets
      bestSubset <- which.min(sds) 
      # this gives the LTS=MCD subset
      idx   <- x_order[bestSubset:(bestSubset + h - 1)] 
      # ranks of LTS subset
      mu    <- means[bestSubset] # LTS=MCD location
      sigma <- sds[bestSubset]   # LTS=MCD scale
      result[i] <- -0.5 * h * log(2 * pi) -
        0.5 * h * log(sigma^2) -
        0.5  / sigma^2 * sum( (xyj[idx] - mu)^2) + 
        (lambda - 1) * sum(sign(x[idx])*log(1 + abs(x[idx])))
    } # ends loop over lambdagrid
    xt <- YJ(x, lambda = lambdagrid[which.max(result)])$yt
    zt <- (xt - median(xt)) / mad(xt)
  } else { # Box-Cox
    for (i in seq_len(length(lambdagrid))) {
      lambda  <- lambdagrid[i]
      xyj     <- BC(x, lambda)$yt
      nbSubsets <- n - h + 1
      sds   <- rep(1, nbSubsets)
      means <- rep(0, nbSubsets)
      for (j in seq_len(nbSubsets)) { # loop over subsets
        datatemp <- xyj[x_order[j:(j + h - 1)]]
        # subset of ordered xyj
        means[j] <- mean(datatemp)
        sds[j]   <- sd(datatemp)
      } # ends loop over subsets
      bestSubset <- which.min(sds) 
      # this gives the LTS=MCD subset
      idx   <- x_order[bestSubset:(bestSubset + h - 1)] 
      # ranks of LTS subset
      mu    <- means[bestSubset] # LTS=MCD location
      sigma <- sds[bestSubset]   # LTS=MCD scale
      result[i] <- -n / 2 * log(sigma^2) + (lambda - 1) * sum(log(x))
    }
    xt <- BC(x, lambdagrid[which.max(result)])$yt
    zt <- (xt - median(xt)) / mad(xt)
  }
  return(list(lambda = lambdagrid[which.max(result)],
              lambdagrid = lambdagrid, LL = result,
              xt = xt, zt = zt))
}


localFivenum = function(x, na.rm = TRUE) {
  # Our version of stats::fivenum where now the quartiles
  # are order statistics, with symmetric ranks.
  # Assumes x is _sorted_ already.
  xna <- is.na(x)
  if (any(xna)) {
    if (na.rm) 
      x <- x[!xna]
    else return(rep.int(NA, 5))
  }
  n <- length(x)
  if (n == 0) 
    rep.int(NA, 5)
  else {
    med <- if ((n %% 2) == 0) {
      (x[(n/2)] + x[(n/2) + 1]) / 2
    } else {
      x[(n - 1)/2 + 1]
    }
    result <- c(x[1], x[ceiling(n / 4)], med,
                x[n - ceiling(n / 4) + 1], x[n])
    result
  }
}


BCr = function(y, lambda, chg = NULL, stdToo = TRUE) {
  # Rectified Box Cox transformation.
  # Like the Classical Box Cox transformation but with a linear 
  # tail on the shrinking side, and a continuous derivative.
  #
  # Arguments:
  # y      : univariate data
  # lambda : transformation parameter for BC
  # chg    : change point: where the linear part begins
  #          when lambda < 1, or ends when lambda > 1.
  #          If NULL, uses the getChangepoint function
  # stdToo : also output poststandardized transformed data.
  #
  if(min(y) <= 0) stop("Data values should be strictly positive")
  yt = rep(NA,times = length(y))
  chgt <- NULL
  if(lambda == 1){ # is already linear
    yt = y-1
    chgt = 0
  }
  if (lambda > 1){
    if(is.null(chg)) { chg = getChangepointBCr(y, lambda = lambda) }
    indl = which(y < chg)
    if (length(indl) > 0) {
      yt[-indl] = (y[-indl]^lambda-1)/lambda
      chgt      = (chg^lambda-1)/lambda
      yt[indl]  = chgt + (y[indl]-chg)*chg^(lambda-1)
    } else {
      yt = (y^lambda-1)/lambda
    }
  } else if(lambda < 1){
    if(is.null(chg)) { chg = getChangepointBCr(y,  lambda = lambda) }
    indu = which(y > chg)    
    if (lambda == 0) {
      if (length(indu) > 0) {
        yt[-indu] = log(y[-indu])
        chgt      = log(chg)
        yt[indu]  = chgt + (y[indu]-chg)/chg
      } else {
        yt = log(y)
      }
    } else {
      if (length(indu) > 0) {
        yt[-indu] = (y[-indu]^lambda-1)/lambda
        chgt      = (chg^lambda-1)/lambda
        yt[indu]  = chgt + (y[indu]-chg)*chg^(lambda-1)
      } else {
        yt = (y^lambda-1)/lambda
      }
    }
  }
  if (stdToo) {
    if (length(y) > 1) {
      locScale <- estLocScale(matrix(yt, ncol = 1), 
                              type = "hubhub")
      zt <- (yt - locScale$loc)/locScale$scale
    } else {
      zt <- yt
    }
  } else {
    zt = NULL
  }
  return(list(yt=yt, chg=chg, chgt=chgt, zt=zt))
}


YJr = function(y, lambda, chg = NULL, prec = 1e-10, stdToo = TRUE) {
  # Rectified Yeo-Johnson transformation.
  # Like the classical YJ transformation but with a linear tail 
  # on the shrinking side, and a continuous derivative.
  #
  # Arguments:
  # y      : univariate data. Assumed shifted to median zero.
  # lambda : transformation parameter for YJ.
  # chg    : change point: where the linear part begins
  #          when lambda < 1, or ends when lambda > 1.
  #          If NULL, uses Q3 for lambda > 1 and
  #          Q1 for lambda > 1.
  # stdToo : also output poststandardized transformed data.
  #
  quarts = localFivenum(y)[2:4]
  yt = rep(NA,times = length(y))
  if (lambda == 1){ # is already linear
    yt = y
    chgt = 0
  }
  indneg = which(y < 0)
  indpos = which(y >= 0)
  #
  if(lambda > 1){
    if(is.null(chg)) { chg = getChangepointYJr(y, lambda = lambda) }
    indl = which(y < chg) # to be linearized
    indb = which(y < 0 & y >= chg) # in between
    yt[indpos] = ((1 + y[indpos])^(lambda) - 1)/lambda
    #
    if(lambda == 2){
      yt[indb] = -log(1-y[indb])
      chgt     = -log(1-chg)
      yt[indl] = chgt + (y[indl]-chg) * YJprimex(chg,lambda)
    }
    else {
      yt[indb] = -((1-y[indb])^(2-lambda)-1)/(2-lambda)
      chgt     = -((1-chg)^(2-lambda)-1)/(2-lambda)
      yt[indl] = chgt + (y[indl]-chg) * YJprimex(chg,lambda)
    }
  }  
  if(lambda < 1){
    if(is.null(chg)) { chg = getChangepointYJr(y, lambda = lambda) }
    # if(chg <= 0) stop("chg should be positive")
    indu = which(y > chg) # to be linearized
    indb = which(y >= 0 & y <= chg) # in between
    yt[indneg] = -((1-y[indneg])^(2-lambda)-1)/(2-lambda)
    #
    if(lambda == 0){
      yt[indb] = log(1 + y[indb])
      chgt     = log(1 + chg)
      yt[indu] = chgt + (y[indu]-chg)*YJprimex(chg,lambda)
    }
    else {
      yt[indb] = ((1+y[indb])^(lambda)-1)/lambda
      chgt     = ((1+chg)^(lambda)-1)/lambda
      yt[indu] = chgt + (y[indu]-chg)*YJprimex(chg,lambda)
    }
  } 
  if (stdToo) {
    if (length(y) > 1) {
      locScale <- estLocScale(matrix(yt, ncol = 1), 
                              type = "hubhub")
      zt <- (yt - locScale$loc) / locScale$scale
    } else { 
      zt <- yt
    } 
  } else {
    zt <- NULL
  }
  return(list(yt=yt, chg=chg, chgt=chgt, zt=zt))
}


rhoBiweight <- function(x, c = 0.5) {
  xbw <- c * (1 - (1 - (x/c)^2)^3)
  xbw[which(abs(x) > c)] <- c
  return(xbw)
}


robnormality <- function(y, b = 0.5) {
  # Assess normality of a dataset after robustly standardizing.
  # Outliers only have a limited effect on this criterion.
  # y is assumed to be _sorted_.
  #
  # args: 
  #   y : vector of observations
  #   b : tuning parameter of biweight rho
  # 
  sy       <- y[!is.na(y)]
  n        <- length(sy)
  locScale <- estLocScale(matrix(sy, ncol = 1), 
                          type = "hubhub")
  sy       <- (sy - locScale$loc) / locScale$scale
  # Theoretical quantiles and differences:
  hprobs     <- ((seq_len(n)) - 1/3)/(n + 1/3)
  theoQs     <- qnorm(hprobs)
  diffs      <- sy - theoQs
  # apply Biweight rho
  rhodiffs <- rhoBiweight(diffs, b)
  crit     <- mean(abs(rhodiffs))
  return(crit)
}