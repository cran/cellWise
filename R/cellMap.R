cellMap <- function(R, indcells = NULL, indrows = NULL, 
                    standOD = NULL, showVals = NULL, D = NULL, 
                    rowlabels = NULL, columnlabels = NULL, 
                    mTitle = "cell map", rowtitle = "cases", 
                    columntitle = "variables", 
                    showrows = NULL, showcolumns = NULL, 
                    nrowsinblock = NULL, 
                    ncolumnsinblock = NULL, 
                    manualrowblocksizes = NULL,
                    manualcolumnblocksizes = NULL,
                    rowblocklabels = NULL,
                    columnblocklabels = NULL,
                    autolabel = TRUE, # but will be ignored
                    columnangle = 90, sizemain = 2, 
                    sizetitles = 1.1, adjustrowlabels = 1,
                    adjustcolumnlabels = 1, colContrast = 1,
                    outlyingGrad = TRUE, 
                    darkestColor = sqrt(qchisq(0.999, 1)), 
                    drawCircles = FALSE) {
  
  # Explanation of the added and modified arguments:
  #
  # rowlabels : Labels of the rows of the matrix R. If NULL, 
  #             these labels are taken as rownames(R), and 
  #             failing that they are 1:nrow(R).
  # columnlabels : Labels of the columns of the matrix R. If NULL,
  #                these labels are taken as colnames(R), and 
  #                failing that they are 1:ncol(R).
  # nrowsinblock : How many rows are combined in a block. Defaults 
  #                to NULL, which asks not to block rows. The 
  #                argument nrowsinblock is overruled by the 
  #                argument manualrowblocksizes when the latter
  #                is specified.
  # ncolumnsinblock : How many columns are combined in a block.
  #                   Defaults to NULL, which asks not to block 
  #                   columns. The argument ncolumnsinblock is
  #                   overruled by the argument 
  #                   manualcolumnblocksizes when the latter is
  #                   specified.
  # manualrowblocksizes : This allows the user to specify their 
  #     own row blocks, unlike the argument nrowsinblock which 
  #     makes all row blocks the same length. The argument takes 
  #     the form c(a,b,...) where a is the length of the first 
  #     block, b is the length of the second, and so on. The 
  #     numbers a,b,... must be strictly positive integers, adding 
  #     up to at most nrow(R). They cannot all be 1, which would 
  #     mean no blocking of rows. Defaults to NULL.
  # manualcolumnblocksizes : Analogous to manualrowblocksizes but 
  #     for columns. It is allowed for one of them to be NULL while
  #     the other is not.
  # rowblocklabels : This allows the user to specify labels for the 
  #     row blocks, whether obtained from nrowsinblock or from
  #     manualrowblocksizes. Defaults to NULL. Will throw an error
  #     if the number of row labels does not match the number of 
  #     blocks.
  # columnblocklabels : Analogous to rowblocklabels but for columns. 
  #     It is allowed for one of them to be NULL while the other
  #     is not.
  
  colorBlocks = function(Xin, rowblocksizes, columnblocksizes, 
                         colContrast) {
    # Auxiliary function, which mixes the colors of the cells
    # in each block.
    n = length(rowblocksizes)
    d = length(columnblocksizes)
    Xblock = matrix(0, nrow = n, ncol = d)
    Xblockgrad = matrix(0, nrow = n, ncol = d)
    rind = cumsum(c(0,rowblocksizes)); rind
    cind = cumsum(c(0,columnblocksizes)); cind
    for (i in seq_len(n)) {
      for (j in seq_len(d)) {
        # select the cells of Xin in the block (i,j):
        Xsel = Xin[(rind[i]+1):rind[i+1],
                   (cind[j]+1):cind[j+1]]
        seltable = tabulate(Xsel, nbins = 4)
        if (sum(seltable) > 0) {
          indmax = which(seltable == max(seltable))[1]
          cntmax = seltable[indmax]
          ncells = rowblocksizes[i] * columnblocksizes[j] 
          gradmax = (cntmax/ncells)^(1/colContrast)
        }
        else {
          indmax = 0
          gradmax = 1
        }
        Xblock[i, j] = indmax
        Xblockgrad[i, j] = gradmax
      }
    }
    return(list(X = Xblock, Xgrad = Xblockgrad))
  }
  
  isVectorInteger = function(x, tol = .Machine$double.eps^0.5){
    y = abs(x - round(x)) < tol 
    sum(1-y) == 0
  }
  
  variable <- rownr <- rescaleoffset <- x <- y <- NULL
  type = "cell" # type="residual" replaces yellow by white
  n = nrow(R)
  d = ncol(R)
  if (!is.null(showVals)) {
    if (!showVals %in% c("D", "R")) {
      stop(paste("Invalid \"showVals\" argument. Should be one of: NULL, \"D\", \"R\""))
    }
  }
  if (is.null(showVals)) {
    D = R
  } else if (showVals == "D" & is.null(D)) {
    stop("When showVals=\"D\" you must input the argument D")
  }
  if (is.null(D)) {
    D = R
  } else { # D exists
    if (!all(dim(D) == dim(R))) 
      stop("The dimensions of D and R must match")
  }
  if (is.null(indcells)) {
    indcells = which(abs(R) > sqrt(qchisq(0.99, 1)))
  }
  
  if(!is.null(rowlabels)){ # rowlabels are given
    if (length(rowlabels) != n) {
      stop(paste0("Number of rowlabels does not match n = ",n))
    }
  } else { # no rowlabels are given
    if (is.null(rownames(R))) {
      rowlabels = seq_len(n)
    } else { rowlabels = rownames(R) }  
  } 
  if(!is.null(columnlabels)){ # columnlabels are given
    if (length(columnlabels) != d) {
      stop(paste0("Number of columnlabels does not match d = ",d))
    }
  } else { # no columnlabels are given
    if (is.null(colnames(R))) {
      columnlabels = seq_len(d)
    } else { columnlabels = colnames(R) }
  }
  # From here on we have n rowlabels and d columnnames. 
  # They are actually case names and variable names, but 
  # we do not rename them which would confuse existing
  # users.
  
  if (!is.null(showcolumns) | !is.null(showrows)) {
    if (is.null(showrows)) {
      showrows = seq_len(n)
    }
    else {
      if (!(all(showrows %in% seq_len(n))))
        stop(" showrows goes out of bounds")
    }
    if (is.null(showcolumns)) {
      showcolumns = seq_len(d)
    }
    else {
      if (!(all(showcolumns %in% seq_len(d))))
        stop(" showcolumns goes out of bounds")
    }
    tempMat = matrix(0, n, d)
    tempMat[indcells] = 1
    tempMat = tempMat[showrows, showcolumns]
    indcells = which(tempMat == 1)
    tempVec = rep(0, n)
    tempVec[indrows] = 1
    tempVec = tempVec[showrows]
    indrows = which(tempVec == 1)
    rm(tempMat, tempVec)
    R = R[showrows, showcolumns]
    D = D[showrows, showcolumns]
    rowlabels = rowlabels[showrows]
    columnlabels = columnlabels[showcolumns]
    n = nrow(R)
    d = ncol(R)
    if (!is.null(standOD))
      standOD = standOD[showrows]
  }
  
  blockRows = blockColumns = FALSE
  if(!is.null(manualrowblocksizes)) {
    if(!is.null(nrowsinblock)){
      cat(paste0("Input argument manualrowblocksizes ",
                 "has overruled argument nrowsinblock.\n"))
      blockRows = TRUE # the rows will be blocked
    }
    msg = paste0(
      "manualrowblocksizes should be a vector with strictly \n",
      "  positive integers, adding up to at most n = ",n)
    if(!is.vector(manualrowblocksizes)) stop(msg)
    if(!is.numeric(manualrowblocksizes)) stop(msg)
    if(!isVectorInteger(manualrowblocksizes)) stop(msg)
    if(sum(manualrowblocksizes < 1) > 0) stop(msg)
    if(sum(manualrowblocksizes) > n) stop(msg)
    if(sum(manualrowblocksizes != 1) == 0 ) {
      stop("All manualrowblocksizes are 1") }
    blockRows = TRUE # the rows will be blocked
  } else { # manualrowblocksizes was not specified
    if(!is.null(nrowsinblock)){
      if(nrowsinblock > 1){
        if (nrowsinblock > n) stop(paste0(
          "Input argument nrowsinblock cannot be ",
          "more than n = ",n))
        blockRows = TRUE # the rows will be blocked
      }
    }
  }
  
  if(!is.null(manualcolumnblocksizes)) {
    if(!is.null(ncolumnsinblock)){
      cat(paste0("Input argument manualcolumnblocksizes ",
                 "has overruled argument ncolumnsinblock.\n"))
    }
    msg = paste0(
      "manualcolumnblocksizes should be a vector with strictly \n",
      "  positive integers, adding up to at most d = ",d)
    if(!is.vector(manualcolumnblocksizes)) stop(msg)
    if(!is.numeric(manualcolumnblocksizes)) stop(msg)
    if(!isVectorInteger(manualcolumnblocksizes)) stop(msg)
    if(sum(manualcolumnblocksizes < 1) > 0) stop(msg)
    if(sum(manualcolumnblocksizes) > d) stop(msg)
    if(sum(manualcolumnblocksizes != 1) == 0 ) {
      stop("All manualcolumnblocksizes are 1") }
    blockColumns = TRUE # the columns will be blocked
  } else {
    if (!is.null(ncolumnsinblock)){
      if(ncolumnsinblock > 1){
        if(ncolumnsinblock > d) stop(paste0(
          "Input argument ncolumnsinblock cannot be ",
          "more than d = ",d)) 
        blockColumns = TRUE # the columns will be blocked
      }
    } 
  }
  
  if((blockRows | blockColumns) & !is.null(showVals)) {
    warning(paste0(
      "The option showVals=\"D\" or showVals=\"R\" cannot be\n", 
      "combined with blocking rows and/or columns,\n", 
      "so showVals is set to NULL here."))
    showVals = NULL
  }
  
  if (type == "residual") outlyingGrad = 1
  X = matrix(0, n, d)
  Xrow = matrix(0, n, 1) # one entry per row, for circles
  Xrow[indrows, 1] = 3   # on the right hand side
  if (type == "cell" | blockRows | blockColumns) { 
    pcells = indcells[indcells %in% which(R >= 0)]
    ncells = indcells[indcells %in% which(R < 0)]
  } else {
    pcells = which(R >= 0)
    ncells = which(R < 0)
  }
  X[ncells] = 1
  X[pcells] = 2
  X[is.na(R)] = 4
  
  if (blockRows | blockColumns) { # starts blocked situation
    rowblocksizes = rep(1,n) # NULL
    if(blockRows){
      if(!is.null(manualrowblocksizes)){
        rowblocksizes = manualrowblocksizes
        n = length(rowblocksizes)
      } else {
        if(!is.null(nrowsinblock)){
          if(nrowsinblock > 1){
            n = floor(n/nrowsinblock)
            rowblocksizes = rep(nrowsinblock,n)
          }  
        }
      }
    }
    columnblocksizes = rep(1,d) # NULL 
    if (blockColumns) {
      if(!is.null(manualcolumnblocksizes)){
        columnblocksizes = manualcolumnblocksizes
        d = length(columnblocksizes)
      } else {
        if(!is.null(ncolumnsinblock)){
          if(ncolumnsinblock > 1){
            d = floor(d/ncolumnsinblock)
            columnblocksizes = rep(ncolumnsinblock,d)
          }
        }
      }
    }
    
    # For coloring the blocks inside the matrix:
    result = colorBlocks(X, rowblocksizes, columnblocksizes,
                         colContrast)
    X = result$X
    Xgrad = result$Xgrad
    if(drawCircles){
      # For coloring (in grey) the circles on the right:
      result = colorBlocks(Xrow, rowblocksizes, c(1), colContrast)
      Xrowgrad = result$Xgrad
      Xrowgrad[result$X == 0] = 0
    }
    
    # Now label the blocks as if they were cells:
    if (blockRows){
      if(is.null(rowblocklabels)){
        cat(paste0("No rowblocklabels were given, so they ", 
                   "are constructed automatically.\n"))
        laby = rowlabels
        rowlabels = rep(0, n) # this is the new (smaller) n
        rind = cumsum(c(0,rowblocksizes))
        for (i in seq_len(n)) {
          if(rowblocksizes[i] == 1){
            rowlabels[i] = laby[rind[i]+1]
          } else {
            rowlabels[i] = paste0(
              laby[rind[i]+1],"-",laby[rind[i+1]])
          }
        }
      } else { # the user has given rowblocklabels
        if (length(rowblocklabels) != n) {
          stop(paste0("The number of rowblocklabels is ", 
                      length(rowblocklabels), 
                      " but there are ",n," row blocks."))
        }        
        rowlabels = rowblocklabels
      }
    }
    if (blockColumns) {
      if(is.null(columnblocklabels)){
        cat(paste0("No columnblocklabels were given, so they ",
                   "are constructed automatically.\n"))
        labx = columnlabels
        columnlabels = rep(0, d)
        cind = cumsum(c(0,columnblocksizes)); cind
        for (j in seq_len(d)) {
          if(columnblocksizes[j] == 1){
            columnlabels[j] = labx[cind[j]+1]
          } else {
            columnlabels[j] = paste0(
              labx[cind[j]+1],"-",labx[cind[j+1]])
          }
        }
      } else { # the user has given columnblocklabels
        if (length(columnblocklabels) != d) {
          stop(paste0("The number of columnblocklabels is ", 
                      length(columnblocklabels), 
                      " but there are ",d," column blocks."))
        }
        columnlabels = columnblocklabels
      }        
    } 
    # From here on the blocks are treated as cells,
    # with lower n and d, and with their own labels.
    
    Xdf = data.frame(cbind(seq(1, n, 1), X))
    colnames(Xdf) = c("rownr", seq(1, d, 1))
    rownames(Xdf) = NULL
    Xdf$rownr = with(Xdf, reorder(rownr, seq(n, 1, -1)))
    mX = reshape2::melt(Xdf, id.var = "rownr", value.name = "CatNr")
    Xgraddf = data.frame(cbind(seq(1, n, 1), Xgrad))
    colnames(Xgraddf) = c("rownr", seq(1, d, 1))
    rownames(Xgraddf) = NULL
    Xgraddf$rownr = with(Xgraddf, reorder(rownr, seq(n, 1, 
                                                     -1)))
    mXgrad = reshape2::melt(Xgraddf, id.var = "rownr", value.name = "grad")
    mX$grad = mXgrad$grad
    mX$rescaleoffset = mXgrad$grad + 10 * mX$CatNr
    if(drawCircles){
      mXrow = data.frame(rownr = seq_len(n), 
                         rescaleoffset = Xrowgrad + 10 * 3)
    }
    scalerange = c(0, 1)
    gradientends = scalerange + rep(c(0, 10, 20, 30, 40), 
                                    each = 2)
    if (type == "cell") 
      colorends = c("yellow", "yellow", "yellow", 
                    "blue", "yellow", "red", "white", 
                    "black", "yellow", "white")
    if (type == "residual") 
      colorends = c("white", "white", "white", 
                    "blue", "white", "red", "white", 
                    "black", "white", "white")
    # ends blocked situation
  } else { # no blocking
    Ddf = data.frame(cbind(seq(1, n, 1), D))
    colnames(Ddf) = c("rownr", seq(1, d, 1))
    rownames(Ddf) = NULL
    Ddf$rownr = with(Ddf, reorder(rownr, seq(n, 1, -1)))
    mD = reshape2::melt(Ddf, id.var = "rownr")
    Rdf = data.frame(cbind(seq(1, n, 1), R))
    colnames(Rdf) = c("rownr", seq(1, d, 1))
    rownames(Rdf) = NULL
    Rdf$rownr = with(Rdf, reorder(rownr, seq(n, 1, -1)))
    mR = reshape2::melt(Rdf, id.var = "rownr")
    Xdf = data.frame(cbind(seq(1, n, 1), X))
    colnames(Xdf) = c("rownr", seq(1, d, 1))
    rownames(Xdf) = NULL
    Xdf$rownr = with(Xdf, reorder(rownr, seq(n, 1, -1)))
    mX = reshape2::melt(Xdf, id.var = "rownr", value.name = "CatNr")
    if (!is.null(showVals)) {
      if (showVals == "D") 
        mX$data = mD$value
      if (showVals == "R") 
        mX$data = mR$value
    }
    if (!outlyingGrad) {
      mX$rescaleoffset = 10 * mX$CatNr
      scalerange = c(0, 1)
      gradientends = scalerange + rep(c(0, 10, 20, 30, 
                                        40), each = 2)
      gradientends
      colorends = c("yellow", "yellow", "blue", 
                    "blue", "red", "red", "white", 
                    "black", "white", "white")
    } else { # if outlyingGrad
      Xgrad = matrix(NA, n, d)
      if (type == "cell") {
        Xgrad[indcells] = abs(R[indcells])
        limL = sqrt(qchisq(0.9, 1))
      }
      else {
        Xgrad = abs(R)
        limL = 0
      }
      limH = darkestColor
      Xgrad[Xgrad > limH] = limH
      Xgrad = ((Xgrad - limL)/(limH - limL))^colContrast
      Xgrad[is.na(Xgrad)] = 0
      Xgraddf = data.frame(cbind(seq(1, n, 1), Xgrad))
      colnames(Xgraddf) = c("rownr", seq(1, d, 1))
      rownames(Xgraddf) = NULL
      Xgraddf$rownr = with(Xgraddf, reorder(rownr, seq(n, 
                                                       1, -1)))
      mXgrad = reshape2::melt(Xgraddf, id.var = "rownr", value.name = "grad")
      mX$grad = mXgrad$grad
      mX$rescaleoffset = mXgrad$grad + 10 * mX$CatNr
      scalerange = c(0, 1)
      gradientends = scalerange + rep(c(0, 10, 20, 30, 
                                        40), each = 2)
      if (type == "cell") 
        colorends = c("yellow", "yellow", 
                      "yellow", "blue", "yellow", 
                      "red", "white", "black", 
                      "white", "white")
      if (type == "residual") 
        colorends = c("white", "white", "white", 
                      "blue", "white", "red", "white", 
                      "black", "white", "white")
    } # ends outlyingGrad
    
    if(drawCircles){
      tempVec = rep(0, n)
      tempVec[indrows] = 1
      mXrow = data.frame(rownr = seq_len(n), rescaleoffset = 40 - 
                           (10 * tempVec))
      rm(tempVec)
      if (is.null(standOD)) {
        mXrow$rescaleoffset[indrows] = mXrow$rescaleoffset[indrows] + 
          1
      } else { # if there is a standOD
        limL = 1
        limH = 3
        standOD[standOD > limH] = limH
        standOD = ((standOD - limL)/(limH - limL))^colContrast
        mXrow$rescaleoffset[indrows] = mXrow$rescaleoffset[indrows] + 
          standOD[indrows]
      } # ends standOD
    } # ends drawCircles
  } # ends unblocked situation
  
  if (drawCircles) { # Xrow is only used here
    circleFun = function(centerx, centery, r, npoints) {
      tt = seq(0, 2 * pi, length.out = npoints)
      xx = centerx + r * cos(tt)
      yy = centery + r * sin(tt)
      return(c(xx, yy))
    }
    columnlabels = c(columnlabels, "", "") 
    # else it ends in the plot with NA NA when d+2
    # or still NA now it was changed to d+1
    centerx = d + 1
    centery = n:1
    radius = 0.4
    npoints = 100
    circlePoints = mapply(circleFun, centerx, centery, radius, 
                          npoints)
    positions = data.frame(rownr = rep(seq_len(n), each = npoints), 
                           x = c(circlePoints[seq_len(npoints), ]), y = c(circlePoints[(npoints + 
                                                                                          1):(2 * npoints), ]))
    datapoly = merge(mXrow, positions, by = c("rownr"))
  } # ends if(drawCircles)
  
  rowlabels = rev(rowlabels)
  base_size = 10
  
  ##### Here ggplot starts
  
  ggp = ggplot(data = mX, aes(variable, rownr)) + {
    geom_tile(aes(fill = rescale(rescaleoffset, from = range(gradientends))), 
              color = "white")  } +
    { if (drawCircles){ 
      geom_polygon(data = datapoly, aes(x = x, y = y, fill = rescale(rescaleoffset, from = range(gradientends)), group = rownr), colour = "black")} } + 
    scale_fill_gradientn(colours = colorends, values = rescale(gradientends), 
                         rescaler = function(x, ...) x, oob = scales::squish) + coord_fixed() + theme_classic(base_size = base_size * 
                                                                                                                1) + labs(x = columntitle, y = rowtitle) + 
    { if(drawCircles){ scale_x_discrete(expand = c(0,0), limits = as.factor(seq(1, d+1 , 1)), labels = columnlabels) } else { scale_x_discrete(expand = c(0,0), limits = as.factor(seq(1, d , 1)), labels = columnlabels) } } + 
    scale_y_discrete(expand = c(0, 0), labels = rowlabels) + 
    ggtitle(mTitle) + 
    theme(legend.position = "none", axis.ticks = element_blank(), 
          plot.title = element_text(size = base_size * sizemain,
                                    hjust = 0.5, vjust = 1, face = "bold"),
          axis.text.x = element_text(size = base_size * 1.8, 
                                     angle = columnangle, hjust = adjustcolumnlabels, 
                                     vjust = 0.5, colour = "black"), axis.text.y = element_text(size = base_size * 
                                                                                                  1.8, angle = 0, hjust = adjustrowlabels, colour = "black"), 
          axis.title.x = element_text(colour = "black", 
                                      size = base_size * sizetitles, vjust = 1), axis.title.y = element_text(colour = "black", 
                                                                                                             size = base_size * sizetitles, vjust = 0), axis.line.x = element_blank(), 
          panel.border = element_blank()) + annotate(geom = "segment", 
                                                     x = 0.5, xend = d + 0.5, y = 0.5, yend = 0.5) + annotate(geom = "segment", 
                                                                                                              x = 0.5, xend = d + 0.5, y = n + 0.5, yend = n + 0.5) + 
    annotate(geom = "segment", x = d + 0.5, xend = d + 
               0.5, y = 0.5, yend = n + 0.5) 
  if (!is.null(showVals)) {
    txtcol = mX$CatNr
    txtcol[txtcol == 0] = "black"
    txtcol[txtcol == 1] = "white"
    txtcol[txtcol == 2] = "white"
    if (type == "residual") {
      txtcol[] = "black"
      txtcol[mXgrad$grad > 0.5] = "white"
    }
    txtcol[txtcol == 4] = "black"
    ggp = ggp + geom_text(aes(label = ifelse(is.na(data), 
                                             sprintf("%1.0f", data), round(data, 1))), size = base_size * 
                            0.5, colour = txtcol, na.rm = TRUE)
  }
  return(ggp)
}