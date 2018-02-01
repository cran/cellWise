
cellMap <- function(D, R, indcells, indrows, showVals = NULL, xlabels = "",
                   ylabels = "", mTitle = "", xtitle = "", ytitle = "",
                   xshowindex = NULL, yshowindex = NULL, xblocksize = 1,
                   yblocksize = 1, autolabel = TRUE, anglex = 90, sizexy = 1.1,
                   hjustXlabels = 1, hjustYlabels = 1, colContrast = 1,
                   outlyingGrad = FALSE) {
  # Draws a cellmap, possibly of a subset of rows and columns of the data,
  # and possibly combining cells into blocks. 
  # The inputs are:
  #
  # D            : the data matrix (required input argument)
  # R            : matrix of cell residuals (required input argument)
  # indcells     : indices of outlying cells (required input argument)
  # indrows      : indices of outlying rows (required input argument)
  # showVals     : show entries ("D"= data matrix, "R"= cell residuals)
  # xlabels      : labels for the x-axis
  # ylabels      : labels for the y-axis
  # mTitle       : main title of the cellMap  
  # xtitle       : title for the x-axis
  # ytitle       : title for the y-axis
  # xshowindex   : indices of the cells that will be shown, in the x direction
  # yshowindex   : indices of the cells that will be shown, in the y direction  
  # xblocksize   : size of combination blocks in the x direction
  # yblocksize   : size of combination blocks in the y direction
  # autolabel    : automatically combines labels of cells in blocks.  
  #                If FALSE, you must provide the final xlabels and/or ylabels.
  # anglex       : angle of the labels on the x-axis
  # sizexy       : size of title for x-axis and y-axis
  # hjustXlabels : adjust x-labels: 0=left, 0.5=centered, 1=right
  # hjustYlabels : adjust y-labels: 0=left, 0.5=centered, 1=right
  # colContrast  : increase contrast with value from 1-5 (1= default) 
  # outlyingGrad : adjust colors in function of outlyingness (only for blocksizes=1)
  
  funcSqueeze <- function(Xin, n, d, xblocksize, yblocksize, colContrast) {
    # function to combine cells into blocks
    Xblock     <- matrix(0, nrow = n, ncol = d)
    Xblockgrad <- matrix(0, nrow = n, ncol = d)
    for (i in 1:n) {
      for (j in 1:d) {
        Xsel     <- Xin[(1 + ((i - 1) * yblocksize)):(i * yblocksize),
                        (1 + ((j - 1) * xblocksize)):(j * xblocksize)]
        seltable <- tabulate(Xsel, nbins = 3)
        cnt0     <- (yblocksize * xblocksize) - sum(seltable)
        if (sum(seltable) > 0) {
          indmax  <- which(seltable == max(seltable))[1]
          cntmax  <- seltable[indmax]
          gradmax <- (cntmax / (xblocksize * yblocksize)) ^ (1 / colContrast)
        } else {
          indmax  <- 0
          gradmax <- 1
        }
        Xblock[i, j]     <- indmax
        Xblockgrad[i, j] <- gradmax
      }
    }
    return(list(X = Xblock, Xgrad = Xblockgrad))
  }
  
  n <- nrow(R)
  d <- ncol(R)
  
  #check input arguments 
  blockMap <- FALSE
  if (xblocksize > 1 | yblocksize > 1 ) {
    blockMap <- TRUE
    if (xblocksize > d) stop('Input argument xblocksize cannot be larger than d')
    if (yblocksize > n) stop('Input argument yblocksize cannot be larger than n')
    if (!is.null(showVals)) warning('The option showVals=D or showVals=R cannot be combined
                                    with xblocksize or yblocksize greater than 1,
                                    so showVals is set to NULL here.')
    showVals <- NULL
  }
  
  if (!blockMap) {
    if (!all(dim(R) == dim(D))) stop('Dimensions of D and R must match')
  }
  
  if (!(blockMap & autolabel == FALSE)) {
    if (length(xlabels) > 0 & length(xlabels) != d) {
      stop(paste('Number of xlabels does not match d = ', d, sep = ""))}
    if (length(ylabels) > 0 & length(ylabels) != n) {
      stop(paste('Number of ylabels does not match n = ', n, sep = ""))}
  }
  
  if (!is.null(showVals)) {
    if (!showVals %in% c("D", "R")) {
      stop(paste("Invalid \"showVals\" argument. Should be one of: NULL, \"D\", \"R\""))
    }
  }
  
  if (!outlyingGrad %in% c(FALSE, TRUE)) {
    stop(paste("Invalid \"outlyingGrad\" argument. Should be TRUE or FALSE"))
  }
  
  if (!(is.null(xshowindex) & is.null(yshowindex))) {
    # here we extract the rows and columns that will be shown.
    if (is.null(xshowindex)) { xshowindex <- 1:d } else {
      if (!(all(xshowindex %in% 1:d))) stop(" xshowindex goes out of bounds")}
    if (is.null(yshowindex)) { yshowindex <- 1:n } else {
      if (!(all(yshowindex %in% 1:n))) stop(" yshowindex goes out of bounds")}
    
    tempMat <- matrix(0, n, d)
    tempMat[indcells] <- 1
    tempMat <- tempMat[yshowindex, xshowindex]
    indcells <- which(tempMat == 1)
    
    tempVec <- rep(0, n)
    tempVec[indrows] <- 1
    tempVec <- tempVec[yshowindex]
    indrows <- which(tempVec == 1) # also works if indrows was empty
    rm(tempMat, tempVec)
    
    if (!blockMap) {
      D <- D[yshowindex, xshowindex] 
    }
    R <- R[yshowindex, xshowindex]
    xlabels <- xlabels[xshowindex]
    ylabels <- ylabels[yshowindex]
    n <- nrow(R)
    d <- ncol(R)
  }  
  
  # create the matrix which indicates the color of each cell
  # 0=yellow, 1=blue, 2=red, 3=black, 4=white
  X <- matrix(0, n, d)
  X[indrows, ] <- 3
  pcells <- indcells[indcells %in% which(R >= 0)]
  ncells <- indcells[indcells %in% which(R < 0)]
  X[ncells] <- 1
  X[pcells] <- 2
  if (!blockMap) X[is.na(D)] <- 4
  
  if (blockMap) { # in this case the input D will be ignored  
    n <- floor(n / yblocksize)
    d <- floor(d / xblocksize)
    
    # create new X and Xgrad
    result <- funcSqueeze(X, n, d, xblocksize, yblocksize, colContrast)
    X      <- result$X
    Xgrad  <- result$Xgrad
    
    # if autolabel=F, labels{x,y} will be used for the blocks.
    if (autolabel == TRUE) { # automatically combine labels for blocks
      
      if (xblocksize > 1 & length(xlabels) > 0) {
        labx <- xlabels
        xlabels <- rep(0, d)
        for (ind in 1:d) {
          xlabels[ind] <- paste(labx[(1 + ((ind - 1) * xblocksize))], "-",
                                labx[(ind * xblocksize)], sep = "")
        }
      }
      
      if (yblocksize > 1 & length(ylabels) > 0) {
        laby <- ylabels
        ylabels <- rep(0, n)
        for (ind in 1:n) {
          ylabels[ind] <- paste(laby[(1 + ((ind - 1) * yblocksize))], "-",
                                laby[(ind * yblocksize) ])
        }
      }
    } else {
      if (length(xlabels) > 0 & length(xlabels) != d) {
        stop(paste(' autolabel=FALSE and number of xlabels is ',
                   length(xlabels), ' but should be ', d, sep = ""))}
      if (length(ylabels) > 0 & length(ylabels) != n) {
        stop(paste(' autolabel=FALSE and number of ylabels is ',
                   length(ylabels), ' but should be ', n, sep = ""))}      
    }
    
    # Melt data matrices for cellMap
    Xdf <- data.frame(cbind(seq(1, n, 1), X))
    colnames(Xdf) <- c("rownr", seq(1, d, 1))
    rownames(Xdf) <- NULL
    Xdf$rownr <- with(Xdf, reorder(rownr, seq(n, 1, -1)))
    mX <- melt(Xdf, id.var = "rownr", value.name = "CatNr") 
    
    Xgraddf <- data.frame(cbind(seq(1, n, 1), Xgrad))
    colnames(Xgraddf) <- c("rownr", seq(1, d, 1))
    rownames(Xgraddf) <- NULL
    Xgraddf$rownr <- with(Xgraddf, reorder(rownr, seq(n, 1, -1)))
    mXgrad <- melt(Xgraddf, id.var = "rownr", value.name = "grad") 
    
    # Combine melted data
    mX$grad <- mXgrad$grad
    mX$rescaleoffset <- mXgrad$grad + 10 * mX$CatNr
    scalerange <- c(0,1)
    gradientends <- scalerange + rep(c(0, 10, 20, 30), each = 2)
    colorends <- c("yellow", "yellow", "yellow", "blue", "yellow",
                   "red", "yellow", "black")
    
  } else {# no blockMap
    
    # Melt data matrices for cellMap
    Ddf <- data.frame(cbind(seq(1, n, 1), D))
    colnames(Ddf) <- c("rownr", colnames(D))
    rownames(Ddf) <- NULL
    Ddf$rownr <- with(Ddf, reorder(rownr, seq(n, 1, -1)))
    mD <- melt(Ddf, id.var = "rownr") 
    
    Rdf <- data.frame(cbind(seq(1, n, 1), R))
    colnames(Rdf) <- c("rownr", colnames(D))
    rownames(Rdf) <- NULL
    Rdf$rownr <- with(Rdf, reorder(rownr, seq(n, 1, -1)))
    mR <- melt(Rdf, id.var = "rownr") 
    
    Xdf <- data.frame(cbind(seq(1, n, 1), X))
    colnames(Xdf) <- c("rownr", seq(1, d, 1))
    rownames(Xdf) <- NULL
    Xdf$rownr <- with(Xdf, reorder(rownr, seq(n, 1, -1)))
    mX <- melt(Xdf, id.var = "rownr", value.name = "CatNr") 
    
    if (!is.null(showVals)) {
      # Combine melted data
      if (showVals == "D") mX$data <- mD$value
      if (showVals == "R") mX$data <- mR$value
    }
    
    if (!outlyingGrad) {
      mX$rescaleoffset <- 10 * mX$CatNr
      scalerange <- c(0, 1)
      gradientends <- rep(c(0, 10, 20, 30, 40), each = 2)
      colorends <- c("yellow", "yellow", "blue", "blue", "red", "red",
                     "black", "black", "white", "white")  
    } else {
      Xgrad <- abs(R)
      Xgrad[Xgrad > 4] <- 4
      Xgrad <- (Xgrad / 4) ^ colContrast
      Xgrad[is.na(Xgrad)] <- 0
      
      Xgraddf <- data.frame(cbind(seq(1, n, 1), Xgrad))
      colnames(Xgraddf) <- c("rownr", seq(1, d, 1))
      rownames(Xgraddf) <- NULL
      Xgraddf$rownr <- with(Xgraddf, reorder(rownr, seq(n, 1, -1)))
      mXgrad <- melt(Xgraddf, id.var = "rownr", value.name = "grad") 
      
      mX$grad <- mXgrad$grad
      mX$rescaleoffset <- mXgrad$grad + 10 * mX$CatNr
      scalerange <- c(0, 1)
      gradientends <- scalerange + rep(c(0, 10, 20, 30, 40), each = 2)
      colorends <- c("yellow", "yellow", "yellow", "blue", "yellow", "red",
                     "black", "black", "white", "white")  
    }
  }
  
  ylabels <- rev(ylabels) # will be plotted from bottom to top
  base_size <- 10
  variable <- rownr <- rescaleoffset <- NULL
  
  ggp <- ggplot(data = mX, aes(variable, rownr))
  
  if (blockMap) {
    ggp <- ggp + geom_tile(aes(fill = 
                                 rescale(rescaleoffset,
                                         from = range(c(0, 1) +
                                                        rep(c(0, 10, 20, 30), 
                                                            each = 2)))),
                           colour = "white")
  }
  if (!blockMap & outlyingGrad) { 
    ggp <- ggp + geom_tile(aes(fill = 
                                 rescale(rescaleoffset,
                                         from = range(c(0,1) +
                                                        rep(c(0, 10, 20, 30, 40), 
                                                            each = 2)))),
                           colour = "white")
  }
  if (!blockMap & !outlyingGrad) {
    ggp <- ggp + geom_tile(aes(fill = rescale(rescaleoffset,
                                              from = c(0, 40))),
                           colour = "white")
  }
  
  ggp <- ggp + scale_fill_gradientn(colours = colorends,
                                    values = rescale(gradientends),
                                    rescaler = function(x, ...) x,
                                    oob = identity) +  
    ggtitle(mTitle) +
    coord_fixed() + 
    theme_classic(base_size = base_size*1) + 
    labs(x = xtitle, y = ytitle) + 
    scale_x_discrete(expand = c(0, 0), labels = xlabels) +
    scale_y_discrete(expand = c(0, 0), labels = ylabels) +
    theme(legend.position = "none",axis.ticks = element_blank(), 
          plot.title = element_text(size = base_size * 2, hjust = 0.5,
                                    vjust = 1, face = "bold"),
          axis.text.x = element_text(size = base_size * 1.8, angle = anglex,
                                     hjust = hjustXlabels, vjust = 0.5, 
                                     colour = "black"),
          axis.text.y = element_text(size = base_size * 1.8, angle = 0,
                                     hjust = hjustYlabels, colour = "black"),
          axis.title.x = element_text(colour = "black",
                                      size = base_size * sizexy, vjust = 1),
          axis.title.y = element_text(colour = "black",
                                      size = base_size * sizexy, vjust = 0)) 
  
  if (!is.null(showVals)) {
    txtcol <- mX$CatNr
    txtcol[txtcol == 0] <- "black"
    txtcol[txtcol == 1] <- "white"
    txtcol[txtcol == 2] <- "white"
    txtcol[txtcol == 3] <- "white"
    txtcol[txtcol == 4] <- "black"
    # the following line places `NA' in missing cells
    # the size of the symbols is also specified here
    ggp <- ggp + geom_text(aes(label = ifelse(is.na(data), sprintf("%1.0f", data),
                                              round(data, 1))),
                           size = base_size * 0.5, colour = txtcol, na.rm = TRUE) 
  } 
  
  return(ggp)
}

