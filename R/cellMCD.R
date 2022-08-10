
cellMCD <- function(X, alpha = 0.75, quant = 0.99, crit = 1e-4, 
                   noCits = 100, lmin = 1e-4,
                   checkPars = list()) {
  #
  # checkPars is as in cellWise::transfo(). If not coreOnly,
  # we run checkDataSet() like we did in other functions.
  #
  #
  #
  # Arguments:
  # X         : A n by d data matrix or data frame to be 
  #             analyzed. Its columns are the variables.
  # alpha     : In each column, at least n*alpha cells must
  #             remain unflagged. Defaults to 75%, should not
  #             be set (much) lower.
  # quant     : Determines the cutoff value to flag cells. 
  #             Defaults to 0.99
  # crit      : The iteration stops when successive covariance
  #             matrices (of the standardized data) differ by
  #             less than crit. Defaults to 1e-4.
  # noCits    : The maximal number of C-steps used.
  # lmin      : a lower bound on the eigenvalues of the 
  #             estimated covariance matrix on the 
  #             standardized data. Defaults to 1e-04.
  #             Should not be smaller than 1e-6.
  # lmax      : if not NULL, an upper bound on the eigenvalues 
  #             of the estimated covariance matrix on the 
  #             standardized data. Could e.g. be set to
  #             2*ncol(X) since the largest eigenvalue is at
  #             most d max_{ij} |C_{ij}| = d max_j |C_{jj}|.
  # checkPars : Optional list of parameters used in the call 
  #            to cellMCD. The options are:
  #     - coreOnly: If TRUE, skip the execution of checkDataset. 
  #         Defaults to FALSE
  #     - numDiscrete:  A column that takes on numDiscrete or 
  #         fewer values will be considered discrete and not
  #         retained in the cleaned data. Defaults to 5.
  #     - precScale: Only consider columns whose scale is larger 
  #         than precScale. Here scale is measured by the median 
  #         absolute deviation. Defaults to 1e-12.
  #     - silent: Whether or not the function progress messages 
  #         should be printed. Defaults to FALSE.
  
  # First some auxiliary functions that are only used here, so 
  # they don't show up in the list of functions:
  
  lmax <- NULL # upper bound on maximum eigenvalue. not used
  
  Cstep <- function(X, W, mu, Sigma, Sigmai, lambdas, h, noits, 
                   precScale=1e-12, nWupdates=1){
    # Performs the C-step for cellMCD.
    # Part 1 updates W for fixed mu and Sigma,
    # Part 2 updates mu and Sigma for fixed W.
    #
    dims <- dim(X)
    n    <- dims[1]
    d    <- dims[2]
    #
    # Part 1: start by updating W. 
    for (jjj in seq_len(nWupdates)) {
      W <- updateW_cpp(X = X, W = W, mu = mu,
                       Sigma = Sigma, Sigmai = Sigmai,
                       lambda = lambdas, h = h)
    }
    
    # Part 2: now execute one EM-step. 
    #
    # First take one E-step.
    # This can be done faster: we should only need to invert 
    # a single matrix here!
    # We could also loop over patterns instead of rows.
    #
    Ximp <- X # will contain suff stats for estimating mu
    bias <- matrix(0, d, d) # more suff stats for estimating Sigma
    #
    for (i in 1:n) { # loops over rows, not patterns
      mis <- which(W[i, ] == 0)
      obs <- which(W[i, ] == 1)
      x <- X[i, ]
      if (length(mis) > 0)  {
        if (length(mis) == d) { # whole row is missing
          ximp <- mu # just plug in old mu
          bias <- bias + Sigma # "bias" is updated by old Sigma
        } else{
          ximp <- x
          Sigmai_temp <- solve(Sigmai[mis, mis]) 
          ximp[mis] <- mu[mis] - Sigmai_temp %*%
            Sigmai[mis, obs, drop = F] %*%
            t(X[i, obs, drop = F] - mu[obs]) 
          bias[mis, mis] <- bias[mis, mis] + Sigmai_temp
        }
        Ximp[i, ] <- ximp
      }
    }
    #
    # Now one M-step, with the same formulas as if the sufficient 
    # statistics Ximp and "bias" came from complete data:
    #
    mu    <- colMeans(Ximp)
    bias  <- bias / n
    Sigma <- cov(Ximp) * (n - 1) / n + bias
    return(list(W = W, mu = mu, Sigma = Sigma))
  }
  
  iterMCD <- function(X,
                     initEst,
                     alpha=0.75,
                     lambdas,
                     precScale=1e-12,
                     crit=1e-4,
                     noCits=100,
                     noEMits=1,
                     nWupdates=1,
                     lmin = NULL,
                     lmax = NULL,
                     silent) {
    #
    if (!is.list(initEst)) stop("initEst should be a list") 
    h <- round(alpha * dim(X)[1])
    p <- ncol(X)
    n <- nrow(X)
    
    mu     <- initEst$mu
    Sigma  <- initEst$Sigma
    Sigmai <- mpinv(Sigma)$Inv
    if (is.null(initEst$W)) { W <- matrix(1, n, p)
    } else { W <- initEst$W }
    W[is.na(X)] <- 0 # to deal with NA's in data
    
    convcrit <- 1
    nosteps  <- 0
    objvals  <- rep(NA, noCits + 1)
    penalty  <- sum(lambdas * colSums(1 - W))
    objvals[nosteps + 1] <- Objective_cpp(X, W, mu, Sigma, Sigmai) + penalty
    if (!silent) cat(paste0("\nObjective at step ", nosteps, " = ",
                           round(objvals[nosteps + 1], 5), "\n"))
    prevW     <- W
    prevmu    <- mu
    prevSigma <- Sigma
    
    while (convcrit > crit && nosteps < noCits) {
      Cresult  <- Cstep(X = X, W = W, mu = mu, Sigma = Sigma,
                       Sigmai = Sigmai, lambdas = lambdas, h = h, 
                       noits = noEMits, nWupdates = nWupdates)
      convcrit <- max(abs((Cresult$Sigma - Sigma))) 
      W       <- Cresult$W
      mu      <- Cresult$mu
      Sigma   <- Cresult$Sigma
      Sigma   <- truncEig(Sigma, lmin, lmax)
      Sigmai  <- solve(Sigma)
      penalty <- sum(lambdas * colSums(1 - W))
      objvl   <- Objective_cpp(X, W, mu, Sigma, Sigmai) + penalty 
      if (objvl > objvals[nosteps + 1]) { 
        return(list(W = prevW, mu = prevmu, Sigma = prevSigma,
                    nosteps = nosteps))
      }
      if (!silent) cat(paste0("Objective at step ", nosteps + 1, " = ",
                              round(objvl, 5), "\n"))
      nosteps <- nosteps + 1
      objvals[nosteps+1] <- objvl
      prevW     <- W
      prevmu    <- mu
      prevSigma <- Sigma
    }
    return(list(W = W, mu = mu, Sigma = Sigma, nosteps = nosteps))
  }
  
  # Here the main function starts:
  #
  X <- as.matrix(X) 
  if (!"coreOnly" %in% names(checkPars)) {
    checkPars$coreOnly <- FALSE
  }
  if (!"numDiscrete" %in% names(checkPars)) {
    checkPars$numDiscrete <- 5
  }
  if (!"precScale" %in% names(checkPars)) {
    checkPars$precScale <- 1e-12
  }
  if (!"silent" %in% names(checkPars)) {
    checkPars$silent <- FALSE
  }
  
  if (!"fracNA" %in% names(checkPars)) {
    checkPars$fracNA <- 0.5
  }
  
  if (!checkPars$coreOnly) {
    X <- checkDataSet(X, fracNA = 0.5, 
                      numDiscrete = checkPars$numDiscrete, 
                      precScale = checkPars$precScale, 
                      silent = checkPars$silent)$remX
  }
  #
  # Robustly standardize the data
  #
  locsca  <- cellWise::estLocScale(X)
  rscales <- locsca$scale
  Xs      <- scale(X, center = locsca$loc, scale = rscales)
  #
  # Check whether there are not too many bad cells:
  #
  Z        <- Xs
  cutf     <- sqrt(qchisq(quant, df = 1))
  margfrac <- colSums(abs(Z) > cutf, na.rm = TRUE) / nrow(Z)
  if (max(margfrac) > 1 - alpha) {
    cat(paste0("\nAt least one variable of X has more than ",
               "100*(1-alpha)% = ", 100 * (1 - alpha), "%",
               "\nof marginal outliers.",
               "\nThe percentages per variable are:\n"))
    print(round(100 * margfrac, 2))
    stop("Too many marginal outliers.")
  }
  Z[which(abs(Z) > cutf)] <- NA
  badfrac <- colMeans(is.na(Z)) # also includes real NAs in X
  if (max(badfrac) > 1 - alpha) {
    cat(paste0("\nAt least one variable of X has more than ",
               "100*(1-alpha)% = ", 100 * (1 - alpha), "%",
               "\nof marginal outliers plus NA's.",
               "\nThe percentages per variable are:\n"))
    print(round(100 * badfrac, 2))
    stop("Too many marginal outliers plus NA's.")
  }
  #
  if (nrow(X) < 5 * ncol(X)) {
    warning(paste0("There are fewer than 5 cases per dimension",
                   " in this data set.\n",
                   "It is not recommended to run cellMCD",
                   " on these data.\n",
                   "Consider reducing the number of variables."))
  }
  if (lmin < 1e-6) stop("lmin should be at least 1e-6.")
  #
  DDCWout <- DDCWcov(Xs, maxCol = 1 - alpha, lmin = lmin, lmax = lmax)
  initEst <- list(mu = DDCWout$center, Sigma = DDCWout$cov)
  #
  # We set the marginal outliers to NA so they stay flagged:
  Xs[abs(Xs) > 3] = NA
  #
  logC    <- -log(diag(solve(initEst$Sigma)))
  cutoff  <- qchisq(quant, df = 1)
  lambdas <- cutoff + logC + log(2 * pi)
  temp    <- iterMCD(Xs, initEst = initEst, alpha = alpha, lambdas = lambdas, 
                     precScale = checkPars$precScale, crit = crit, 
                     noCits = noCits, noEMits = 1, nWupdates = 1, 
                     lmin = lmin, lmax = lmax, 
                     silent = checkPars$silent)
  W <- temp$W
  rownames(W) <- rownames(X)
  colnames(W) <- colnames(X)
  out   <- allpreds_cpp(Xs, temp$Sigma, temp$mu, W)
  preds <- out$preds
  if (sum(is.na(as.vector(preds))) > 0 ) stop("There are missing preds!")
  cvars <- out$cvars
  if (sum(is.na(as.vector(cvars))) > 0 ) stop("There are missing cvars!") 
  if (min(as.vector(cvars)) <= 0) stop("There are cvars <= 0 !")
  rownames(preds) <- rownames(cvars) <- rownames(X)
  colnames(preds) <- colnames(cvars) <- colnames(X)
  S <- diag(rscales) %*% temp$Sigma %*% diag(rscales)
  if (!checkPars$silent) {
    percflag <- 100 * colMeans(1 - W, na.rm = TRUE)
    cat("Percentage of flagged cells per variable:\n")
    print(round(percflag,2))
  }
  mu    <- locsca$loc + temp$mu * rscales
  preds <- scale(preds, center = FALSE, scale = 1 / rscales)
  preds <- scale(preds, center = -locsca$loc, scale = FALSE)
  csds  <- scale(sqrt(cvars), center = FALSE, scale = 1 / rscales)
  Ximp  <- X
  Ximp[which(W == 0)] <- preds[which(W == 0)]
  Zres  <- (X - preds) / csds
  locsc <- cellWise::estLocScale(Zres, center = FALSE)
  Zres  <- scale(Zres, center = FALSE, scale = locsc$scale)
  return(list(mu = mu, S = S,
              W = W, preds = preds,
              csds = csds, Ximp = Ximp,
              Zres = Zres, rscales = rscales,
              nosteps = temp$nosteps,
              X = X, quant = quant))
} 


plot_cellMCD = function(cellout, type = "Zres/X", whichvar = NULL,
                        horizvar = NULL, vertivar = NULL,  
                        hband = NULL, vband = NULL, drawellipse = T,
                        opacity = 0.5, identify = FALSE, 
                        ids=NULL, labelpoints = T, vlines = FALSE,
                        clines = TRUE, main = NULL,
                        xlab = NULL, ylab = NULL, xlim = NULL,
                        ylim = NULL, cex = 1, cex.main = 1.2, 
                        cex.txt = 0.8, cex.lab = 1, line=2.0){
  #
  # Function for making plots based on cellMCD output.
  #
  # Arguments:
  # cellout      output of function cellMCD()
  # type         "index", "Zres/X", "Zres/pred", "X/pred", or
  #              "bivariate".
  # whichvar     number or name of the variable to be plotted.  
  #              Not applicable when type == "bivariate".
  # horizvar     number or name of the variable to be plotted on the 
  #              horizontal axis. Only when type == "bivariate".
  # vertivar     number or name of the variable to be plotted on the 
  #              vertical axis. Only when type == "bivariate".
  # hband        draw a horizontal tolerance band? TRUE or FALSE. 
  #              NULL yields TRUE for types "index", "Zres/X", 
  #              and "Zres/pred".
  # vband        draw a vertical tolerance band? TRUE or FALSE.
  #              NULL yields TRUE for types "Zres/X", "Zres/pred", 
  #              and "X/pred".
  # drawellipse  whether to draw a 99% tolerance ellipse. Only
  #              for type == "bivariate".
  # opacity      opacity of the plotted points: 1 is fully opaque,
  #              less is more transparent.
  # identify     if TRUE, identify cases by mouseclick, then Esc.
  # ids          vector of case numbers to be emphasized in the plot.
  #              If NULL or of length zero, none are emphasized.
  # labelpoints  if TRUE, labels the points in ids by their
  #              row name in X.
  # vlines       for the points in ids, draw dashed vertical lines 
  #              from their standardized residual to 0 when type is
  #              "index", "Zres/X", or "Zres/pred". Draws dashed
  #              vertical ines to the diagonal for type "X/pred". 
  #              Can be TRUE or FALSE, default is FALSE.
  # clines       only for type == "bivariate". If TRUE, draws
  #              a red connecting line from each point in ids to 
  #              its imputed point, shown in blue.
  #
  # The following arguments are plot options, to finalize plots
  # for presentation:
  #
  # main         main title of the plot. If NULL, it is constructed
  #              automatically from the arguments.  
  # xlab         overriding label for x-axis, unless NULL.
  # ylab         overriding label for y-axis, unless NULL.
  # xlim         overriding limits of horizontal axis.
  # ylim         overriding limits of vertical axis.
  # cex          size of plotted points.
  # cex.main     size of the main title.
  # cex.lab      size of the axis labels.
  # cex.txt      size of the point labels.
  # line         distance of axis labels to their axis.
  #
  # Invisible output:
  # out      NULL, except when identify == TRUE. Then a list with:
  #          $ids    : the case number(s) that were identified
  #          $coords : coordinates of all points in the plot.
  
  # First some auxiliary functions:
  #
  identfy = function(xcoord,ycoord){
    # identify points in a plot(x,y)
    coordinates <- cbind(xcoord,ycoord)
    message("Press the escape key to stop identifying.")
    iout <- identify(coordinates, order = TRUE)
    ids <- iout$ind[order(iout$order)]
    if (length(ids) > 0) {
      cat("Identified point(s): ")
      print(ids)
    }   
    return(list(ids=ids, coords=coordinates[ids,,drop=F]))
  }
  
  addlabels = function(xcoord, ycoord, ids, labs, cex.txt = 0.8,
                       labtype = NULL){
    # Adds labels to plot. When labs = "i" it is the case number.
    # Also, labs can be the set of row names of the dataset.
    # For labs = "letters" we plot a, b, c, ...  
    #
    len = length(ids)
    if(len == 0) stop("ids has no elements")
    if(is.null(labtype)) { mylabs = labs[ids] 
    } else {
      if(labtype == "i") mylabs = ids
      if(labtype == "letters") mylabs = letters[seq_len(len)]
    }
    text(x = xcoord[ids], y = ycoord[ids], labels = mylabs, 
         cex = cex.txt)
  }
  
  vlines2diag = function(xcoord, ycoord, ids, lty=2, col="red"){
    # For indices in isd, plot vertical residual line 
    for(i in seq_len(length(ids))){ # i=2
      xys = c(xcoord[ids[i]],xcoord[ids[i]],
              ycoord[ids[i]],xcoord[ids[i]])
      lines(matrix(xys, ncol = 2, byrow=F), lty=lty, col=col)
    }
  }
  
  vlines2zero = function(xcoord, ycoord, ids, lty=2, col="red"){
    # For indices in isd, plot vertical residual line 
    for(i in seq_len(length(ids))){ # i=2
      xys = c(xcoord[ids[i]],xcoord[ids[i]],
              ycoord[ids[i]],0)
      lines(matrix(xys, ncol = 2, byrow=F), lty=lty, col=col)
    }
  }
  
  # Replaced ellipse::ellipse() by the simple function below,
  # so we do not have a dependency on library(ellipse).
  ellipsepoints = function(covmat, mu, quant=0.99, npoints = 100)
  { # computes points of the ellipse t(x-mu)%*%covmat%*%(x-mu) = c
    # with c = qchisq(quant,df=2)
    if (!all(dim(covmat) == c(2, 2))) stop("covmat is not 2 by 2")
    eig = eigen(covmat)
    U = eig$vectors
    R = U %*% diag(sqrt(eig$values)) %*% t(U) # square root of covmat
    angles = seq(0, 2*pi, length = npoints+1)
    xy = cbind(cos(angles),sin(angles)) # points on the unit circle
    fac = sqrt(qchisq(quant, df=2))
    scale(fac*xy%*%R, center = -mu, scale=FALSE)
  }  
  
  # Here the main function starts:
  #
  cutf = sqrt(qchisq(cellout$quant, df=1))
  mycol <- adjustcolor("black", alpha.f = opacity)
  redcol <- adjustcolor("red", alpha.f = 1) # = opacity)
  X = as.matrix(cellout$X)
  if(length(dim(X)) != 2) stop("cellout$X is not a matrix")
  ncol = ncol(X)
  nrow = nrow(X)
  cn = colnames(X)
  if(is.null(cn)) cn = seq_len(ncol)
  rn = rownames(X)
  if(is.null(rn)) rn = seq_len(nrow)
  if(type %in% c("index", "Zres/X", "Zres/pred", "X/pred")){
    if(is.null(whichvar)){
      stop("You must specify the variable whichvar to plot.")
    }
    if(whichvar %in% seq_len(ncol)){
      j = whichvar
      varlab = cn[j] 
    } else { 
      if(whichvar %in% cn){
        j = which(cn == whichvar)
        varlab = whichvar
      } else { stop(paste0("whichvar = ",whichvar," is not valid")) }
    }
    flagged = which(cellout$W[,j] == 0) 
    Xj     = X[,j]
    preds  = cellout$preds[,j]
    csds   = cellout$csds[,j]
    Zres   = (Xj - preds)/csds # standardized residuals
    locsc  = cellWise::estLocScale(Zres, center=F)
    Zres   = scale(Zres, center=F, scale=locsc$scale) 
    lspred = cellWise::estLocScale(preds)
    lsX    = cellWise::estLocScale(Xj)
    #
    if(type == "index"){
      ycoord = Zres 
      xcoord = seq_len(length(ycoord))
      hlim = c(-cutf, cutf)
      obsp = which(!is.na(ycoord)) # points to plot
      if(is.null(ylim)) ylim = range(c(ycoord[obsp],hlim), na.rm=T)    
      plot(xcoord[obsp], ycoord[obsp], xlim=xlim, ylim=ylim, pch=16, 
           xlab="", ylab="", col=mycol, cex=cex) 
      if(is.null(xlab)) xlab = "index"
      title(xlab = xlab, line=line, cex.lab = cex.lab)
      if(is.null(ylab)) ylab = paste0("standardized residual of ",varlab)
      title(ylab = ylab, line=line, cex.lab = cex.lab)
      if(is.null(main)) main = 
        paste0("index plot: standardized residual of ",varlab)
      title(main = main, line = 1, cex.main = cex.main)
      if(is.null(hband) || hband==T){
        abline(h = hlim, lwd = 3, col="darkgray") }
      # in an index plot we cannot draw a meaningful vband.
      points(xcoord[flagged], ycoord[flagged], pch=16, col=redcol)
      if(length(ids) > 0){
        if(labelpoints) {
          addlabels(xcoord, ycoord, ids, labs = rn, cex.txt = cex.txt)
        }
        if(vlines == TRUE){
          vlines2zero(xcoord, ycoord, ids, lty=2, col="red")
        }
      }
    }
    #
    if(type == "Zres/X"){
      xcoord = Xj
      ycoord = Zres
      lsx  = cellWise::estLocScale(xcoord) 
      vlim = c(lsx$loc - cutf*lsx$scale, lsx$loc + cutf*lsx$scale)
      hlim = c(-cutf, cutf)
      obsp = which(!is.na(xcoord) & !is.na(ycoord)) # points to plot
      if(is.null(xlim)) xlim = range(c(xcoord[obsp],vlim), na.rm=T)
      if(is.null(ylim)) ylim = range(c(ycoord[obsp],hlim), na.rm=T)    
      plot(xcoord[obsp], ycoord[obsp], xlim=xlim, ylim=ylim, pch=16, 
           xlab="", ylab="", col=mycol, cex=cex) 
      if(is.null(xlab)) xlab=varlab
      title(xlab = xlab, line=line, cex.lab = cex.lab)
      if(is.null(ylab)) ylab = paste0("standardized residual of ",varlab)
      title(ylab = ylab, line=line, cex.lab = cex.lab)
      if(is.null(main)) main = paste0("standardized residual versus X for ",varlab)
      title(main = main,line = 1, cex.main = cex.main) 
      if(is.null(hband) || hband==T){
        abline(h = hlim, lwd = 3, col="darkgray") }  
      if(is.null(vband) || vband==T){
        abline(v = vlim, lwd = 3, col="darkgray") }
      points(xcoord[flagged], ycoord[flagged], pch=16, col=redcol)
      if(length(ids) > 0){
        if(labelpoints) {      
          addlabels(xcoord, ycoord, ids, labs = rn, cex.txt = cex.txt)
        }
        if(vlines == TRUE){
          abline(h=0) 
          vlines2zero(xcoord, ycoord, ids, lty=2, col="red") 
        }
      }
    }
    #
    if(type == "Zres/pred"){  
      xcoord = preds
      ycoord = Zres
      lsx  = cellWise::estLocScale(xcoord) 
      vlim = c(lsx$loc - cutf*lsx$scale, lsx$loc + cutf*lsx$scale)
      hlim = c(-cutf, cutf)
      obsp = which(!is.na(xcoord) & !is.na(ycoord)) # points to plot
      if(is.null(xlim)) xlim = range(c(xcoord[obsp],vlim), na.rm=T)
      if(is.null(ylim)) ylim = range(c(ycoord[obsp],hlim), na.rm=T)    
      plot(xcoord[obsp], ycoord[obsp], xlim=xlim, ylim=ylim, pch=16, 
           xlab="", ylab="", col=mycol, cex=cex) 
      if(is.null(xlab)) xlab = paste0("predicted ",varlab)
      title(xlab=xlab, line=line, cex.lab = cex.lab)
      if(is.null(ylab)) ylab = paste0("standardized residual of ",varlab)
      title(ylab = ylab, line=line, cex.lab = cex.lab)
      if(is.null(main)) main = paste0(
        "standardized residual versus prediction for ",varlab)
      title(main = main,line = 1, cex.main = cex.main)
      if(is.null(hband) || hband==T){
        abline(h = hlim, lwd = 3, col="darkgray") }  
      if(is.null(vband) || vband==T){
        abline(v = vlim, lwd = 3, col="darkgray") }
      points(xcoord[flagged], ycoord[flagged], pch=16, col=redcol)
      if(length(ids) > 0){
        if(labelpoints) {
          addlabels(xcoord, ycoord, ids, labs = rn, cex.txt = cex.txt)
        }
        if(vlines == TRUE){
          abline(h=0) 
          vlines2zero(xcoord, ycoord, ids, lty=2, col="red") 
        }
      }
    }
    #
    if(type == "X/pred"){
      xcoord = preds
      ycoord = Xj
      lsx  = cellWise::estLocScale(xcoord) 
      lsy  = cellWise::estLocScale(ycoord)
      vlim = c(lsx$loc - cutf*lsx$scale, lsx$loc + cutf*lsx$scale)
      hlim = c(lsy$loc - cutf*lsy$scale, lsy$loc + cutf*lsy$scale)
      obsp = which(!is.na(xcoord) & !is.na(ycoord)) # points to plot
      if(is.null(xlim)) xlim = range(c(xcoord[obsp],vlim), na.rm=T)
      if(is.null(ylim)) ylim = range(c(ycoord[obsp],hlim), na.rm=T)    
      plot(xcoord[obsp], ycoord[obsp], xlim=xlim, ylim=ylim, pch=16, 
           xlab="", ylab="", col=mycol, cex=cex) 
      if(is.null(xlab)) xlab = paste0("predicted ",varlab)
      title(xlab = xlab, line=line, cex.lab = cex.lab)
      if(is.null(ylab)) ylab = paste0("observed ",varlab)
      title(ylab = ylab, line=line, cex.lab = cex.lab)
      if(is.null(main)) main = paste0(varlab," versus its prediction")
      title(main = main,line = 1, cex.main = cex.main)
      abline(0,1)
      if(!is.null(hband) && hband == TRUE){
        abline(h = hlim, lwd = 3, col="darkgray") }
      if(!is.null(vband) && vband == TRUE){
        abline(v = vlim, lwd = 3, col="darkgray") }
      points(xcoord[flagged], ycoord[flagged], pch=16, col=redcol)
      if(length(ids) > 0){
        if(labelpoints) {
          addlabels(xcoord, ycoord, ids, labs = rn, cex.txt = cex.txt)
        }
        if(vlines == TRUE){
          vlines2diag(xcoord, ycoord, ids, lty=2, col="red")
        }   
      }
    }
  } else {
    if(type == "bivariate"){
      if(is.null(horizvar)) {
        stop(paste0("You must specify the variable horizvar\n",
                    "to plot on the horizontal axis.")) }
      if(is.null(vertivar)) {
        stop(paste0("You must specify the variable verticar\n",
                    "to plot on the vertical axis.")) }
      if(horizvar %in% seq_len(ncol)){
        jj = horizvar
        hlab = cn[jj] 
      } else { 
        if(horizvar %in% cn){
          jj = which(cn == horizvar)
          hlab = horizvar
        } else { stop(paste0("horizvar = ",horizvar," is not valid")) }
      }
      if(vertivar %in% seq_len(ncol)){
        kk = vertivar
        vlab = cn[kk] 
      } else { 
        if(vertivar %in% cn){
          kk = which(cn == vertivar)
          vlab = vertivar
        } else { stop(paste0("vertivar = ",vertivar," is not valid")) }
      }
      if(jj==kk) stop("horivar and vertivar should differ.")
      xcoord = X[,jj]
      ycoord = X[,kk]
      imp = cellout$Ximp
      ell = ellipsepoints(covmat = cellout$S[c(jj,kk),c(jj,kk)],
                          mu = cellout$mu[c(jj,kk)], quant=0.99,
                          npoints=500)
      lsx = cellWise::estLocScale(xcoord) 
      lsy = cellWise::estLocScale(ycoord)
      vlim = c(lsx$loc - cutf*lsx$scale, lsx$loc + cutf*lsx$scale)
      hlim = c(lsy$loc - cutf*lsy$scale, lsy$loc + cutf*lsy$scale)
      obsp = which(!is.na(xcoord) & !is.na(ycoord)) # points to plot
      if(is.null(xlim)) {
        xlim = range(c(xcoord[obsp],imp[obsp,jj],ell[,1],vlim), na.rm=T)
      }
      if(is.null(ylim)) {
        ylim = range(c(ycoord[obsp],imp[obsp,kk],ell[,2],hlim),na.rm=T)
      }
      flagged = which(cellout$W[,jj]*cellout$W[,kk] == 0) 
      plot(xcoord[obsp], ycoord[obsp], xlim=xlim, ylim=ylim, pch=16, 
           xlab="", ylab="", col=mycol, cex=cex) 
      if(is.null(xlab)) xlab = hlab
      title(xlab = xlab, line=line, cex.lab = cex.lab)
      if(is.null(ylab)) ylab = vlab
      title(ylab = ylab, line=line, cex.lab = cex.lab)
      if(is.null(main)) main = paste0(vlab," versus ",hlab)
      title(main = main,line = 1, cex.main = cex.main)
      if(is.null(hband)) hband = F
      if(hband == T) abline(h = hlim, lwd=3, col="darkgray")
      if(is.null(vband)) vband = F
      if(vband == T) abline(v = vlim, lwd=3, col="darkgray")
      if(drawellipse) lines(ell, lwd=3, col="darkgray")
      points(xcoord[flagged], ycoord[flagged], pch=16, col=redcol)
      if(length(ids) > 0){
        if(labelpoints) {
          addlabels(xcoord, ycoord, ids, labs = rn, cex.txt = cex.txt)
        }
        if(is.null(clines) || clines==T){
          imp = cellout$Ximp
          for(i in ids) {
            if(i %in% obsp){
              lines(x=c(X[i,jj],imp[i,jj]),y=c(X[i,kk],imp[i,kk]), 
                    lty=1, col="red")
              points(x=imp[i,jj],y=imp[i,kk],pch=16,col="blue")
            }
          }
        }   
      }
    } else  stop(paste0("type = \"",type,"\" is invalid"))
  }
  if(identify) invisible(identfy(xcoord,ycoord))
}





truncEig <- function(S, lmin = NULL, lmax = NULL) {
  # Truncates the eigenvalues of S to between lmin and lmax.
  nrS <- nrow(S)
  if (ncol(S) != nrS) stop(" S is not square")
  if (mean(as.vector(abs(S - t(S)))) > 1e-8) {
    stop(" S is not symmetric")
  }
  Sout <- S
  if (!(is.null(lmin) && is.null(lmax))) {
    if (!is.null(lmin)) {
      if (lmin < 0) stop(paste0("lmin = ",lmin," must be >= 0"))
      if (!(lmin < Inf)) stop(paste0("lmin = ",lmin," must be finite"))
    }
    if (!is.null(lmax)) {
      if (lmax < 0) stop(paste0("lmax = ",lmax," must be >= 0"))
      if (!(lmax < Inf)) stop(paste0("lmax = ",lmax," must be finite"))
    }
    eig <- eigen(0.5 * (S + t(S)))
    vals <- eig$values
    if (!is.null(lmin)) { # there is a lower bound
      vals <- pmax(vals, lmin) 
    }
    if (!is.null(lmax)) { # there is an upper bound
      vals <- pmin(vals, lmax) 
    }    
    Sout <- eig$vectors %*% diag(vals) %*% t(eig$vectors) 
  } 
  return(Sout)
}
