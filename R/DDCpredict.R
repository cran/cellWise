

DDCpredict <- function (Xnew, InitialDDC, DDCpars = NULL) 
{ # Only added a check whether the number of columns of Xnew
  # matches that of InitialDDC$remX, to prevent the earlier
  # unintelligible error messages, and removed some useless
  # wnq() lines.
  predictCol2 = function(colj, U, ngbrs, corrweight, robslopes, 
                         combinRule) { # auxiliary function
    j = colj[1]
    colj = colj[-1]
    U = matrix(U[-1, ], ncol = dim(U)[2])
    contributors = (corrweight[j, ] > 0)
    if (length(contributors) < 1) {
      estcol = 0
    }
    else {
      ngb1 = ngbrs[j, contributors]
      slopes1 = robslopes[j, contributors]
      corrwt1 = corrweight[j, contributors]
      ZestAllh = t(t(matrix(U[, ngb1], nrow = dim(U)[1])) * 
                     slopes1)
      if (combinRule == "wmean") {
        estcol = apply(ZestAllh, 1, weighted.mean, w = corrwt1, 
                       na.rm = TRUE)
      }
      if (combinRule == "wmedian") {
        estcol = apply(ZestAllh, 1, weightedMedian, w = corrwt1, 
                       na.rm = TRUE)
      }
      if (combinRule == "mean") {
        estcol = apply(ZestAllh, 1, mean, na.rm = TRUE)
      }
      if (combinRule == "median") {
        estcol = apply(ZestAllh, 1, median, na.rm = TRUE)
      }
      estcol
    }
  }
  
  # Main function starts here
  if (is.null(DDCpars)) {
    DDCpars = InitialDDC$DDCpars
  }
  else {
    if (!is.list(DDCpars)) {
      stop("DDCpars must be a list")
    }
    InitialDDC$DDCpars[names(DDCpars)] <- DDCpars
    DDCpars <- InitialDDC$DDCpars
  }
  if (!"precScale" %in% names(DDCpars)) {
    DDCpars$precScale <- 1e-12
  }
  if (!"tolProb" %in% names(DDCpars)) {
    DDCpars$tolProb <- 0.99
  }
  if (!"corrlim" %in% names(DDCpars)) {
    DDCpars$corrlim <- 0.5
  }
  if (!"combinRule" %in% names(DDCpars)) {
    DDCpars$combinRule <- "wmean"
  }
  if (!"tolProbCell" %in% names(DDCpars)) {
    DDCpars$tolProbCell <- DDCpars$tolProb
  }
  if (!"tolProbRow" %in% names(DDCpars)) {
    DDCpars$tolProbRow <- DDCpars$tolProb
  }
  if (!"includeSelf" %in% names(DDCpars)) {
    DDCpars$includeSelf <- 1
  }
  if (DDCpars$tolProb < 0.5) 
    stop("tolProb must be >= 0.5")
  if (DDCpars$tolProb >= 1) 
    stop("tolProb must be < 1.0")
  qCell = sqrt(qchisq(DDCpars$tolProbCell, 1))
  qRow = sqrt(qchisq(DDCpars$tolProbRow, 1))
  includeSelf = DDCpars$includeSelf
  combinRule = DDCpars$combinRule
  corrlim = DDCpars$corrlim
  uniDetect = limitFilt
  if (is.data.frame(Xnew) | is.matrix(Xnew) | is.vector(Xnew)) {
    Xnew = data.matrix(Xnew)
  }
  else {
    stop("Data matrix must be of class matrix or data.frame")
  }
  if ((nrow(Xnew) > ncol(Xnew)) & ncol(Xnew) == 1) 
    Xnew = t(Xnew)
  n = nrow(Xnew)
  d = ncol(Xnew)
  dold = ncol(InitialDDC$remX) # PR 20221203
  if(d != dold) stop(paste0(
    "Xnew should have ",dold," columns, that correspond\n",
    "  to those of InitialDDC$remX")) # PR 20221203
  locX = InitialDDC$locX
  Z = sweep(Xnew, 2, locX)
  scaleX = InitialDDC$scaleX
  Z = sweep(Z, 2, scaleX, "/")
  indNAs = which(is.na(Z))
  U = uniDetect(Z, qCut = qCell)
  rownames(U) = rownames(Z)
  colnames(U) = colnames(Z)
  UniIndex = setdiff(which(is.na(U)), indNAs)
  nbngbrs = InitialDDC$nbngbrs
  ngbrs = InitialDDC$ngbrs
  robcors = InitialDDC$robcors
  robslopes = InitialDDC$robslopes
  if (includeSelf) 
    corrweight = abs(robcors)
  if (corrlim > 0) {
    corrweight[corrweight < corrlim] = 0
  }
  if (!includeSelf) 
    colStandalone = which(rowSums(corrweight) == 0)
  if (includeSelf) 
    colStandalone = which(rowSums(corrweight[, -1]) == 0)
  colConnected = which(!((seq_len(d)) %in% colStandalone))
  indexStandalone = col2cell(colStandalone, n = n)
  indexStandalone = indexStandalone[indexStandalone %in% UniIndex]
  numiter = 1
  for (iter in seq_len(numiter)) {
    Zest = U
    U = rbind(seq_len(d), U)
    Zest[, colConnected] = apply(U[, colConnected], 2, predictCol2, 
                                 U = U, ngbrs = ngbrs, corrweight = corrweight, robslopes = robslopes, 
                                 combinRule = combinRule)
    U = U[-1, ]
    Zest[, colConnected] = t(InitialDDC$deshrinkage * t(Zest[, 
                                                             colConnected]))
    Zest[is.na(Zest)] = 0
    Zres = Z - Zest
    Zres[, colStandalone] = Z[, colStandalone]
    Zres[, colConnected] = scale(matrix(Zres[, colConnected], 
                                        ncol = length(colConnected)), center = FALSE, InitialDDC$scalestres)
    indcells = which(abs(Zres) > qCell)
    U[indcells] = NA
  }
  rm(U)
  indcells = setdiff(indcells, col2cell(colStandalone, n = n))
  indcells = unique(sort(c(indcells, indexStandalone)))
  Ti = integer(0)
  indrows = integer(0)
  indall = indcells
  compT = function(rowi) {
    mean(pchisq(rowi^2, 1), na.rm = T) - 0.5
  }
  Ti = as.vector(apply(Zres, 1, FUN = compT))
  Ti = scale(Ti, InitialDDC$medTi, InitialDDC$madTi)
  rownames(Ti) = rownames(Z)
  indrows = which(is.na(uniDetect(Ti, qCut = qRow)) & (Ti > 
                                                         0))
  indall = unique(c(indcells, row2cell(indrows, n, d)))
  Zest = sweep(sweep(Zest, 2, scaleX, "*"), 2, locX, 
               "+")
  Xnew[indcells] = Zest[indcells]
  Xnew[indNAs] = Zest[indNAs]
  attr(Zres, "scaled:scale") = NULL
  attr(Ti, "scaled:center") = NULL
  attr(Ti, "scaled:scale") = NULL
  DDCpars <- DDCpars[!(names(DDCpars) %in% c("tolProbCell", 
                                             "tolProbRow", "includeSelf"))]
  return(list(DDCpars = DDCpars, locX = locX, scaleX = scaleX, 
              Z = Z, nbngbrs = nbngbrs, ngbrs = ngbrs, robcors = robcors, 
              robslopes = robslopes, deshrinkage = InitialDDC$deshrinkage, 
              Xest = Zest, scalestres = InitialDDC$scalestres, stdResid = Zres, 
              indcells = indcells, Ti = Ti, medTi = InitialDDC$medTi, 
              madTi = InitialDDC$madTi, indrows = indrows, indNAs = indNAs, 
              indall = indall, Ximp = Xnew))
}

## AUXILIARY FUNCTIONS:


limitFilt <- function(v,qCut) {
  # Detects outliers and sets them to NA.
  # Assumes that the data have already been standardized.
  vout <- v
  vout[(abs(v) > qCut)] <- NA
  return(vout)
}

wnq <- function(string,qwrite=0){ # auxiliary function
  # writes a line without quotes
  if(qwrite==1) write(noquote(string),file="",ncolumns=100)
}


pnq <- function(string,qwrite=0){ # auxiliary function
  # prints a line without quotes
  if(qwrite==1) print(noquote(string))
} 


col2cell = function(colNrs,n) {
  # Transforms column indices to cellwise indices.
  # Here colNrs is a vector with column numbers between 1 and d.
  cindex = t(matrix(rep((colNrs-1)*n,n),ncol=n,byrow=FALSE))
  cindex = cindex + seq(1,n,1) # contains the cell numbers
  return(as.vector(cindex))
}

row2cell = function(rowNrs,n,d) {
  # Transforms row indices to cellwise indices.
  # Here rowNrs is a vector with row numbers between 1 and n.
  as.vector(t(matrix(rep(rowNrs,d),ncol=d,byrow=FALSE))+
              seq(0,n*(d-1),n))
}