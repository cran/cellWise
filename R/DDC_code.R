
DetectDeviatingCells = function(X,DDCpars=list()){ 
  # X is the input data, and must be a matrix or a data frame.
  #   It must always be provided.
  #
  # DDCpars contains all the options chosen by the user, but has
  # a default:
  #
  if(is.null(DDCpars)) { DDCpars=list()}
  if(!is.list(DDCpars)){stop("DDCpars must be a list")}
  if(!"fracNA" %in% names(DDCpars)){DDCpars$fracNA=0.5}
  if(!"numDiscrete" %in% names(DDCpars)){DDCpars$numDiscrete=3}
  if(!"precScale" %in% names(DDCpars)){DDCpars$precScale=1e-12}
  if(!"tolProb" %in% names(DDCpars)){DDCpars$tolProb=0.99}
  if(!"corrlim" %in% names(DDCpars)){DDCpars$corrlim=0.5}
  if(!"combinRule" %in% names(DDCpars)){DDCpars$combinRule="wmean"}
  if(!"includeSelf" %in% names(DDCpars)){DDCpars$includeSelf=TRUE}
  if(!"rowdetect" %in% names(DDCpars)){DDCpars$rowdetect=TRUE}
  if(!"returnBigXimp" %in% names(DDCpars)){DDCpars$returnBigXimp=FALSE}
 
  
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
  fracNA        = DDCpars$fracNA
  numDiscrete   = DDCpars$numDiscrete
  precScale     = DDCpars$precScale
  returnBigXimp = DDCpars$returnBigXimp
  
  # Check the data set and set aside columns and rows that do
  # not satisfy the conditions:
  out1 = checkDataSet(X,fracNA,numDiscrete,precScale)

  # Carry out the actual DetectDeviatingCells algorithm on
  # the remaining dataset out1$remX :
  out = DDCcore(out1$remX,DDCpars)

  if(returnBigXimp){
    Ximp = data.matrix(X,rownames.force = TRUE)
    Ximp[out1$rowInAnalysis,out1$colInAnalysis] = out$Ximp
    out$Ximp = Ximp
  }
  
  return(c(out1,out))
}
  

checkDataSet = function(X,fracNA=0.5,numDiscrete=3,precScale=1e-12){
  # This function checks the dataset X, and sets aside certain
  # columns and rows that do not satisfy the conditions.
  #
  # fracNA      : only consider columns and rows with fewer NAs than this.
  # numDiscrete : a column that takes on numDiscrete or fewer values
  #               will be considered discrete and not used in the analysis.
  # precScale   : only consider columns whose scale is > precScale.
  #               Here scale is measured by the median absolute deviation.
  #
  
  wnq = function(string,qwrite=1){ # auxiliary function
    # writes a line without quotes
    if(qwrite==1) write(noquote(string),file="",ncolumns=100)
  }
  
  pnq = function(string,qwrite=1){ # auxiliary function
    # prints a line without quotes
    if(qwrite==1) print(noquote(string))
  } 
  
  if(!is.data.frame(X) & !is.matrix(X)) {
    stop("The input data must be a matrix or a data frame") }    
  
  n = nrow(X)
  if(n < 3) stop(" The input data must have at least 3 rows (cases)")  
  d = ncol(X)
  if(d < 2) stop(" The input data must have at least 2 columns (variables)")  
  
  wnq(" ")
  wnq(paste(" The input data has ",n, " rows and ",
            d," columns.",sep=""))  
  # Add column names and row names if these were not given:
  if(is.matrix(X)) { X = data.frame(X) } 
  # This keeps the row names and column names if they exist, else creates
  # them as 1, 2, 3, ... for rows and V1, V2, V3,... for columns.
  
  remX = X # remX will be the remaining part of the data
  colInAnalysis = sapply(remX,is.numeric)
  numgoodcol = sum(colInAnalysis)
  vecNotNumeric = (colInAnalysis == FALSE)
  numbadcol = sum(vecNotNumeric) # can be 0
  namesNotNumeric = NULL  
  if(numbadcol > 0) {
    wnq(" ")
    wnq(paste(" The input data contained ",numbadcol,
              " non-numeric columns (variables).",sep=""))
    wnq(" Their column names are:")
    wnq(" ")    
    namesNotNumeric = colnames(remX)[vecNotNumeric]    
    pnq(namesNotNumeric)   
    wnq(" ")    
    if(numgoodcol > 1) {
      wnq(" These columns will be ignored in the analysis.")
      wnq(paste(" We continue with the remaining ",numgoodcol,
                " numeric columns:",sep="")) 
      remX = remX[colInAnalysis]
    } else { 
      if(numgoodcol == 0) stop(" No columns remain, so we stop.")
      if(numgoodcol == 1) stop(" Only 1 column remains, so we stop.")
    }    
    wnq(" ")
    pnq(names(which(colInAnalysis)))
  }
  
  # Turn data into a matrix
  remX = data.matrix(remX,rownames.force = TRUE)
  # This keeps the row names and remaining column names.
  
  # 2. Deselect column(s) containing only the case number
  #
  caseNumber = (1:nrow(remX))
  distFromCaseNumber = function(colj,caseNumber){
    # assumes colj is a vector
    mean(abs(colj-caseNumber),na.rm=F) # can be NA
  }
  dists = apply(remX,2,distFromCaseNumber,caseNumber)
  dists[!is.finite(dists)] = 1 # takes out NA, NaN, Inf, -Inf
  vecbadcol  = (dists == 0)
  numbadcol  = sum(vecbadcol) # can be 0
  goodcol    = (vecbadcol == FALSE)
  numgoodcol = sum(goodcol)
  namesCaseNumber = NULL
  if(numbadcol > 0) {
    wnq(" ")
    wnq(paste(" The data contained ",numbadcol," columns that were",
              " identical to the case number",sep=""))
    wnq(" (number of the row).")
    wnq(" Their column names are:")
    wnq(" ")
    namesCaseNumber = colnames(remX)[vecbadcol]
    pnq(namesCaseNumber)
    wnq(" ")
    if(numgoodcol > 1) {
      wnq(" These columns will be ignored in the analysis.")
      wnq(paste(" We continue with the remaining ",numgoodcol,
                " columns:",sep="")) 
      remX = remX[,goodcol,drop=F]
    } else { 
      if(numgoodcol == 0) stop(" No columns remain, so we stop.")
      if(numgoodcol == 1) stop(" Only 1 column remains, so we stop.")
    }    
    # Update array colInAnalysis using goodcol:
    colInAnalysis[colInAnalysis==TRUE] = goodcol # overwrites correctly
    wnq(" ")
    pnq(names(which(colInAnalysis)))
  }
  
  
  # 3. Deselect variables with over fracNA% of missing values
  #    (e.g. fracNA=0.20). Then update the vector colInAnalysis.
  #
  remX[!is.finite(remX)] = NA # sets NA, NaN, Inf, -Inf all to NA
  acceptNA = nrow(remX)*fracNA
  NAcounts   = colSums(is.na(remX))
  goodcol    = (NAcounts <= acceptNA)
  numgoodcol = sum(goodcol)
  vecNAcol   = (goodcol == FALSE)
  numNAcol   = sum(vecNAcol)
  namesNAcol = NULL
  if(numNAcol > 0) {
    wnq(" ")
    wnq(paste(" The data contained ",numNAcol," columns with over ",
              round(100*fracNA,2),"% of NAs.",sep=""))
    wnq(" Their column names are:")
    wnq(" ")    
    namesNAcol = colnames(remX)[vecNAcol]
    pnq(namesNAcol)    
    wnq(" ")    
    if(numgoodcol > 1) {
      wnq(" These columns will be ignored in the analysis.")
      wnq(paste(" We continue with the remaining ",numgoodcol,
                " columns:",sep=""))
      remX = remX[,goodcol,drop=FALSE]
    } else { 
      if(numgoodcol == 0) stop(" No columns remain, so we stop.")
      if(numgoodcol == 1) stop(" Only 1 column remains, so we stop.")
    }    
    colInAnalysis[colInAnalysis==TRUE] = goodcol
    wnq(" ")
    pnq(names(which(colInAnalysis)))
  }
  
  
  # 4. Deselect rows with too many NAs.
  #    Create the vector rowInAnalysis.
  #
  acceptNA   = ncol(remX)*fracNA
  NAcounts   = rowSums(is.na(remX))
  goodrow    = (NAcounts <= acceptNA)
  numgoodrow = sum(goodrow)
  vecNArow   = (goodrow == FALSE)
  numNArow   = sum(vecNArow)
  rowInAnalysis = goodrow # in case we need to remove more rows later.
  namesNArow = NULL
  if(numNArow > 0) {
    wnq(" ")
    wnq(paste(" The data contained ",numNArow," rows with over ",
              round(100*fracNA,2),"% of NAs.",sep=""))
    wnq(" Their row names are:")
    wnq(" ")
    namesNArow = rownames(remX)[vecNArow]
    pnq(namesNArow)    
    wnq(" ")
    if(numgoodrow > 2) {
      wnq(" These rows will be ignored in the analysis.")
      wnq(paste(" We continue with the remaining ",numgoodrow,
                " rows:",sep="")) 
      remX = remX[goodrow,,drop=F]
    } else { 
      if(numgoodrow == 0) stop(" No rows remain, so we stop.")
      if(numgoodrow == 1) stop(" Only 1 row remains, so we stop.")
      if(numgoodrow == 2) stop(" Only 2 rows remain, so we stop.")      
    }
    wnq(" ")    
    pnq(names(which(rowInAnalysis)))
  }
  
  
  # 5. Deselect discrete variables, loosely defined as variables that
  #    take on numDiscrete or fewer values, such as binary variables.
  #
  countValues = function(colj){
    # assumes colj is a vector
    # length(unique(colj))    # counts NA as a value    
    sum(!is.na(unique(colj))) # only counts non-NAs
  }
  valueCount  = apply(remX,2,countValues)
  goodcol     = (valueCount > numDiscrete)
  numgoodcol  = sum(goodcol)
  vecbadcol   = (goodcol == F)
  numbadcol   = sum(vecbadcol) # can be 0
  namesDiscrete = NULL
  if(numbadcol > 0) {
    wnq(" ")
    wnq(paste(" The data contained ",numbadcol," discrete columns with ",
              numDiscrete," or fewer values.",sep=""))
    wnq(" Their column names are:")
    wnq(" ")
    namesDiscrete = colnames(remX)[vecbadcol]
    pnq(namesDiscrete)    
    wnq(" ")
    if(numgoodcol > 1) {
      wnq(" These columns will be ignored in the analysis.")
      wnq(paste(" We continue with the remaining ",numgoodcol,
                " columns:",sep="")) 
      remX = remX[,goodcol,drop=FALSE]
    } else { 
      if(numgoodcol == 0) stop(" No columns remain, so we stop.")
      if(numgoodcol == 1) stop(" Only 1 column remains, so we stop.")
    }
    colInAnalysis[colInAnalysis==TRUE] = goodcol
    wnq(" ")
    pnq(names(which(colInAnalysis)))
  }
  
  
  # 6. Deselect columns for which the median absolute deviation is
  #    zero. This is equivalent to saying that 50% or more of its
  #    values are equal.
  #
  colScale    = apply(remX,2,mad,na.rm=TRUE)
  goodcol     = (colScale > precScale)
  numgoodcol  = sum(goodcol)
  vecbadcol   = (goodcol == FALSE)
  numbadcol   = sum(vecbadcol) # can be 0
  namesZeroScale = NULL
  if(numbadcol > 0) {
    wnq(" ")
    wnq(paste(" The data contained ",numbadcol," columns with zero",
              " or tiny median absolute deviation.",sep=""))
    wnq(" Their column names are:")
    wnq(" ")
    namesZeroScale = colnames(remX)[vecbadcol]
    pnq(namesZeroScale)
    wnq(" ")
    if(numgoodcol > 1) {
      wnq(" These columns will be ignored in the analysis.")
      wnq(paste(" We continue with the remaining ",numgoodcol,
                " columns:",sep="")) 
      remX = remX[,goodcol,drop=FALSE]
    } else { 
      if(numgoodcol == 0) stop(" No columns remain, so we stop.")
      if(numgoodcol == 1) stop(" Only 1 column remains, so we stop.")
    }    
    colInAnalysis[colInAnalysis==TRUE] = goodcol
    wnq(" ")
    pnq(names(which(colInAnalysis)))
  }
  
  # check whether we have reduced the size of X
  if(nrow(remX) < n | ncol(remX) < d){
    wnq(" ")
    wnq(paste(" The final data set we will analyze has ",
              nrow(remX)," rows and ",ncol(remX),
              " columns.",sep=""))
    wnq(" ")
  }
  
  return(list(colInAnalysis=which(colInAnalysis),
              rowInAnalysis=which(rowInAnalysis),
              namesNotNumeric=namesNotNumeric,
              namesCaseNumber=namesCaseNumber,
              namesNAcol=namesNAcol,
              namesNArow=namesNArow,
              namesDiscrete=namesDiscrete,
              namesZeroScale=namesZeroScale,
              remX=remX))
}


DDCcore = function(X,DDCpars=NULL){
  # This function performs the actual computations of DetectDeviatingCells,
  # and can be called for analyzing real data or in a simulation.
  #
  # X is the input data, and must be a matrix or a data frame.
  # It must always be provided.
  #
  # DDCpars contains all the options chosen by the user, but has a default.
  
  # First, we include all the functions we need here and are not
  # to be used by themselves.
  
  ## FUNCTIONS TO TRANSFORM INDICES IN A MATRIX:
  
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
    as.vector(t(matrix(rep(rowNrs,d),ncol=d,byrow=FALSE))+seq(0,n*(d-1),n))
  }
  
  ## UNIVARIATE ESTIMATORS OF LOCATION AND SCALE:
  
  loc1StepM = function(x,c1=3,precScale=1e-12) {
    # Computes the first step of an algorithm for
    # a location M-estimator using the biweight.
    # Note that c1 is expressed in unnormalized MAD units.
    # In the usual units it is thus c1/qnorm(3/4). 
    x = x[!is.na(x)] # we always take out NAs
    medx = median(x)
    ax = abs(x - medx)
    denom = c1 * median(ax)
    mu = if(!is.na(medx) && denom > precScale) {
      ax = ax/denom
      w = 1 - ax * ax
      w = ((abs(w) + w)/2)^2
      sum(x * w)/sum(w)
    }
    else medx
    return(mu)
  }
  
  scale1StepM = function(x,c2=2.5,delta=0.844472,precScale=1e-12) {
    # Computes the first step of an algorithm for
    # a scale M-estimator using the Huber rho. 
    # Assumes that the univariate data in x has already
    # been centered, e.g. by subtracting loc1StepM.
    # Note that c2 is expressed in unnormalized MAD units.
    # In the usual units it is thus c2/qnorm(3/4).
    # If you change c2 you must also change delta.
    # 
    x = x[!is.na(x)] # we always take out NAs
    n = length(x)
    sigma0 = median(abs(x))
    if(is.na(sigma0) || c2*sigma0 < precScale) { return(sigma0)
    } else {
      x = x/sigma0
      rho = x^2
      rho[rho > c2^2] = c2^2
      return(sigma0 * sqrt(sum(rho)/(n*delta)))
    }
  }
  
  ## UNIVARIATE OUTLIER DETETCTION:
  
  limitFilt = function(v,qCut) {
    # Detects outliers and sets them to NA.
    # Assumes that the data have already been standardized.
    vout = v
    vout[(abs(v) > qCut)] = NA
    return(vout)
  }
  
  equiGYfilt = function(v,qCut,miter = 30) {
    # Detects outliers and sets them to NA. This is a
    # permutation-equivariant version of the GY filter.
    # Assumes that the data have already been standardized.
    qCut = qCut^2
    rawEquiGYfilt = function(v,qCut) {
      # raw version, has no iteration
      v2 = v^2
      n = length(v2)
      u = sort(v2)
      i0 = which(u < qCut) # initial inliers
      n0 = 0
      if (length(i0) > 0) {
        i0 = rev(i0)[1]
        dn = max(pmax(pchisq(u[i0:n],1) - (i0:n - 1)/n, 0))
        n0 = round(dn * n)
      }
      v.na = v
      if (n0 > 0){
        cutoff = u[(n - n0 + 1)]
        v.na[!(v2 < cutoff)] = NA
      } 
      return(as.vector(v.na))
    }
    #
    converge = 0
    iter = 0  
    n = length(v)
    observed = !is.na(v)
    vobs = v[observed]
    nobs = length(vobs)
    id = 1:nobs
    while (converge == 0 & iter < miter) {
      iter = iter + 1
      vobs = rawEquiGYfilt(vobs,qCut=qCut)
      id = id[!is.na(vobs)]
      if (!any(is.na(vobs))) 
        converge = 1
      vobs = na.omit(vobs)
    }
    vobsout = rep(NA,nobs)
    vobsout[id] = vobs
    vout = rep(NA,n)
    vout[observed] = vobsout
    return(vout)
  }
  
  ## FUNCTIONS FOR ROBUST CORRELATION AND SLOPE:
  
  corrGKWLS = function(xcol,qCorr,colj,precScale=1e-12,...) {
    # Computes a robust correlation between the columns xcol and
    # colj using the formula of Gnanadesikan-Kettenring (GK),
    # followed by a Weighted Pearson correlation.
    # qCorr is a quantile used to determine the weights.
    #
    # Assumes colj is a vector with same length as xcol
    # and that all normalizations have already happened.  
    
    if (sum(!is.na(xcol+colj)) <= 3) { #we need at least 3 valids
      corr = 0
    } else {
      corr = (scale1StepM(xcol+colj,precScale=precScale,...)^2 -
                scale1StepM(xcol-colj,precScale=precScale,...)^2)/4
      # This corr should not be NA since data columns with too many NAs
      # have already been taken out. But just in case someone increases
      # the allowed fraction fracNA too much, we add the precaution:  
      if(is.na(corr)) { corr = 0 } else {
        corr = min(0.99,max(-0.99,corr))
        # Compute reweighted Pearson correlation:
        corrMatInv = diag(2) # same as diag(c(1,1))
        corrMatInv[c(2,3)] = -corr # replace off-diagonal elements
        corrMatInv = corrMatInv/abs(1-corr^2) # divide by determinant
        xtemp = cbind(xcol,colj)
        RDi2  = rowSums((xtemp %*% corrMatInv) * xtemp)
        rowSel = (RDi2 < qCorr)
        rowSel[is.na(rowSel)]=FALSE
        if (any(rowSel)) {
          xcol = xcol[rowSel]
          colj = colj[rowSel]
          corr = sum(xcol*colj)/sqrt(sum(xcol^2)*sum(colj^2))  
        } else {
          corr = 0
        }     
      }
    }
    return(corr)
  }
  
  
  kBestCorr = function(colj,U,robCorr=corrGKWLS,k,qCorr,precScale=1e-12){
    # For a given column colj this computes the k highest absolute 
    # correlations (obtained by robCorr) between colj and the columns 
    # of the matrix U, and also reports which k columns are selected.
    #
    # Assumes that the first entry of colj is the number of the column.
    # Assumes the remainder of colj is a vector with same length as the 
    # remainder of the columns of U, and that all normalizations have
    # already happened.
    # Assumes k <= d-1 which is ensured before calling this function.
    j    = colj[1]
    colj = colj[-1] # take out first entry
    U    = U[-1,]   # take out first row
    n    = nrow(U)
    if(!(length(colj) == n)) stop(" colj and U do not match")
    if(k > (ncol(U)-1)) stop(" k is larger than ncol(U)-1")
    allCorrs = apply(U,2,FUN=robCorr,qCorr=qCorr,colj=colj,precScale=precScale) 
    selected = order(abs(allCorrs),decreasing=T)
    notSelf  = !(selected == j)  
    selected = selected[notSelf] # removes correlation of j with j
    # selected has length d-1
    selected = selected[1:k] # always works since k <= d-1
    # We stack indices and correlations into a single vector to make it
    # easier to use apply() later:
    c(selected,allCorrs[selected])
  }
  
  
  slopeMedWLS = function(xcol,qRegr,colj,precScale=1e-12) {
    # Computes the slope of a robust regression without intercept
    # of the column colj on the column xcol.
    # The computation starts by Median regression, followed by
    # weighted least squares (WLS) in which the weights are
    # determined by the quantile qRegr.
    #
    # Assumes that colj is a vector with the same length as xcol
    # and that both columns are already centered.
    
    ratio = colj/xcol
    if (sum(!is.na(ratio)) <= 3) { #we need at least 3 valids
      slope = 0
    } else {
      rawb = median(ratio,na.rm=TRUE) # raw slope
      if(is.na(rawb)) { slope = 0 } else {
        rawb = min(2,max(-2,rawb))
        # Now compute weighted LS slope:
        r = colj - rawb*xcol # raw residuals
        cutoff = qRegr*scale1StepM(r,precScale=precScale) 
        # cutoff can be zero, which is okay.
        rowSel = !(abs(r) > cutoff) # selects the inliers
        rowSel[is.na(rowSel)]=FALSE  
        if (any(rowSel)) {
          yw = colj[rowSel]
          xw = xcol[rowSel] 
          slope = sum(xw*yw)/(sum(xw^2)) # slope of colj ~ xcol + 0
          if(is.na(slope)) { slope = 0 }
        } else {
          slope = 0
        }
      }
    }
    
    return(slope)
  }
  
  
  compSlopes = function(colj,U,robSlope=slopeMedWLS,k,qRegr,precScale=1e-12){
    # For a given column colj this computes the slopes (obtained by
    # robSlope) of regressing colj on each of k given columns of 
    # the matrix U.
    #
    # Assumes that the final k entries of colj are the indices of 
    # those k columns.
    # Assumes the remainder of colj is a vector with same length as the 
    # remainder of the columns of U, and that all of these columns are
    # already centered.
    n    = nrow(U) - k
    ngbr = colj[-(1:n)] # has length k, may contain zeroes
    colj = colj[1:n]
    U    = U[1:n,]
    slopes = rep(0.0,k) # initialize
    ngbr = ngbr[!(ngbr == 0)] # keep only the necessary neighbors
    if(length(ngbr) > 0) {
      if(length(ngbr)==1) {
        b = robSlope(U[,ngbr],qRegr=qRegr,colj=colj,precScale=precScale)   
      } else {
        b = apply(U[,ngbr],2,FUN=robSlope,qRegr=qRegr,colj=colj,
                  precScale=precScale)  
      } 
      slopes[1:length(b)] = b
    }
    return(slopes)
  }
  
  ## ESTIMATING VALUES IN A COLUMN:
  
  predictCol = function(colj,U,ngbrs,corrweight,robslopes,combinRule){
    # Predicts the values in column colj using the set 'ngbrs' of
    # columns of U, by applying the combination rule 'combinRule' whose
    # inputs are the weights in 'corrweight' and the slopes in 'robslopes'.
    #
    # Assumes that the first entry of colj is the number of the column.
    # Assumes the remainder of colj is a vector with same length as the 
    # remainder of the columns of U, and that all of these columns are
    # already centered. 
    j    = colj[1]
    colj = colj[-1] # take out first entry
    U    = U[-1,]   # take out first row
    contributors = (corrweight[j,] > 0)
    if(length(contributors) < 1){
      estcol = 0
    } else {
      ngb1    = ngbrs[j,contributors]
      slopes1 = robslopes[j,contributors]
      corrwt1 = corrweight[j,contributors]
      ZestAllh = t(t(U[,ngb1]) * slopes1) # is n by k matrix
      # Predicts column j from each column h, using slope(Z_j ~ Z_h).
      # This array has the estimates from k variables.    
      if (combinRule =="wmean"){
        estcol = apply(ZestAllh,1,weighted.mean,w=corrwt1,na.rm=TRUE)
      } 
      if (combinRule =="wmedian"){
        estcol = apply(ZestAllh,1,weightedMedian,w=corrwt1,na.rm=TRUE)  
      }
      if (combinRule =="mean"){
        estcol = apply(ZestAllh,1,mean,na.rm=TRUE) 
      }
      if (combinRule =="median"){
        estcol = apply(ZestAllh,1,median,na.rm=TRUE)  
      } 
      estcol
    }    
  }
  
  ## DESHRINKAGE:
  
  deShrink = function(colj,Z,robSlope,qRegr,precScale=1e-12) {
    # Deshrinks the column colj by multiplying it by the robSlope
    # of column j of the matrix Z on colj.
    #
    # Assumes that the first entry of colj is the number of the column.
    # Assumes the remainder of colj is a vector with same length as the 
    # columns of Z, and that both columns were already centered.
    j    = colj[1]
    colj = colj[-1] # the regressor (from Zest)
    zj   = Z[,j]    # column with response (from Z)
    a    = robSlope(colj,qRegr=qRegr,zj,precScale=precScale)
    a * c(0,colj)
  }
  
  wnq = function(string,qwrite=1){ # auxiliary function
    # writes a line without quotes
    if(qwrite==1) write(noquote(string),file="",ncolumns=100)
  }
  
  ## Here the actual computations start:
  # DDCpars contains all the options chosen by the user, but has a default.
  if(is.null(DDCpars)) { 
    DDCpars = list(precScale=1e-12,tolProb=0.99,corrlim=0.5,
                   combinRule ="wmean",includeSelf=T,rowdetect=TRUE)}
  #
  # Retrieve parameters from the list:
  #
  precScale   = DDCpars$precScale
  if(DDCpars$tolProb < 0.5) stop(" tolProb must be >= 0.5")
  if(DDCpars$tolProb >= 1)  stop(" tolProb must be < 1.0")
  qCorr       = qchisq(DDCpars$tolProb,2)
  # = cutoff used for bivariate correlation.
  qRegr       = sqrt(qchisq(DDCpars$tolProb,1))
  # = cutoff used for bivariate slope.
  qCell       = sqrt(qchisq(DDCpars$tolProb,1))
  # = cutoff used for flagging cells.
  qRow        = sqrt(qchisq(DDCpars$tolProb,1)) 
  # = cutoff used for flagging rows.
  includeSelf = DDCpars$includeSelf
  corrlim     = DDCpars$corrlim
  combinRule  = DDCpars$combinRule
  rowdetect   = DDCpars$rowdetect
  
  # Now keeping fixed:
  robLoc     = loc1StepM    # DDCpars$robLoc
  robScale   = scale1StepM  # DDCpars$robScale  
  robSlope   = slopeMedWLS  # DDCpars$robSlope  
  robCorr    = corrGKWLS    # DDCpars$robCorr  
  uniDetect  = limitFilt    # DDCpars$uniDetect
  # robLoc    : robust univariate location estimator.
  # robScale  : robust univariate scale estimator.
  # robSlope  : robust slope estimator.
  # robCorr   : robust correlation estimator.
  # uniDetect : limitFilt compares absolute value to a cutoff.
  #             equiGYfilt means use permutation-equivariant version
  #             of the Gervini-Yohai outlier identifier.
  
  # if(combinRule == "wmedian") { library(matrixStats) }    
  # library matrixStats contains the function weightedMedian
  
  # Turn data into a matrix
  if (is.data.frame(X) | is.matrix(X)) { X = data.matrix(X) } else {
    stop("Data matrix must be of class matrix or data.frame") }
  n = nrow(X)
  d = ncol(X)
  if(d < 2) stop(" DDCcore needs at least 2 columns")

  wnq(paste(" The computation started at: ",date(),sep=""))
  
  #####################################
  ##    STEP 1: STANDARDIZE DATA     ##
  #####################################
  
  # Robust standardization
  locX   = apply(X,2,FUN=robLoc,precScale=precScale)
  Z      = sweep(X,2,locX)
  scaleX = apply(Z,2,FUN=robScale,precScale=precScale)
  Z      = sweep(Z,2,scaleX,"/")
  
  #####################################
  ##    STEP 2: UNIVARIATE ANALYSIS  ##
  #####################################  
  
  missIndex    = which(is.na(X))
  # Univariate analysis by column: replace outliers by NAs
  U = apply(Z,2,uniDetect,qCut=qCell)
  rownames(U) = rownames(Z)
  colnames(U) = colnames(Z)    
  UniIndex    = setdiff(which(is.na(U)),missIndex) #does not include original missings

  #####################################################
  ##    STEP 3: CALCULATE CORRELATIONS AND SLOPES    ##
  #####################################################      
   
  k = min(d-1,1000) # -1 since we do not have to compute the
                    # correlation of a column with itself
 
  # For each column j of U, find the k columns h != j of U that
  # it has the highest absolute correlation robCorr with:
  U  = rbind(1:d,U)
  bc = apply(U,2,kBestCorr,U=U,robCorr=robCorr,k=k,
             qCorr=qCorr,precScale=precScale) 
  U = U[-1,]
  ngbrs   = t(bc[1:k,,drop=F])         # identities of these columns
  robcors = t(bc[(k+1):(2*k),,drop=F]) # the correlations
  
  corrweight = abs(robcors) # should have no NAs
  if (corrlim > 0) { corrweight[corrweight < corrlim] = 0 }
  ngb0 = ngbrs
  ngb0[corrweight == 0] = 0
  U    = rbind(U,t(ngb0))
  robslopes = apply(U,2,compSlopes,U=U,robSlope=robSlope,k=k,
                    qRegr=qRegr,precScale=precScale)
  if(k==1) { robslopes = matrix(robslopes, nrow=1) }
  robslopes = t(robslopes)
  U = U[1:n,] 

  colStandalone = which(rowSums(corrweight)==0)
  colConnected  = which(!((1:d) %in% colStandalone))
  
  indexStandalone = col2cell(colStandalone,n=n)
  indexStandalone = indexStandalone[indexStandalone %in% UniIndex]
  # = list of flagged cells in standalone variables.
  
  if(includeSelf){ 
    # if you want to include column j in its own prediction:
    ngbrs      = cbind(1:d,ngbrs)
    robcors    = cbind(rep(1.0,d),robcors)
    corrweight = cbind(rep(1.0,d),corrweight)
    robslopes  = cbind(rep(1.0,d),robslopes)
  }
  
  numiter = 1 # increase this if you want to iterate
  for (iter in 1:numiter){
  
  ####################################
  ##    STEP 4 : ESTIMATE CELLS     ##
  ####################################     

  Zest = U # These values will remain for standalone columns.  

  # Estimation for connected variables:
  U = rbind(1:d,U)
  Zest[,colConnected] = apply(U[,colConnected],2,predictCol,U=U,ngbrs=ngbrs,
        corrweight=corrweight,robslopes=robslopes,combinRule=combinRule) 
  U = U[-1,]

  ####################################
  ##    STEP 5 : DESHRINKAGE        ##
  ####################################
  
  # Deshrinkage: rescale Zest[,j] using robSlope of Z[,j] on Zest[,j]
  Zest = rbind(1:d,Zest)
  Zest[,colConnected] = apply(Zest[,colConnected],2,deShrink,Z=Z,
               robSlope=robSlope,qRegr=qRegr,precScale=precScale)
  Zest = Zest[-1,]
    
  # Finally, all NAs are replaced by zeroes:
  Zest[is.na(Zest)] = 0
    
  ####################################
  ##    STEP 6 : FLAGGING CELLS     ##
  ####################################      
  
  # Compute cell residuals:
  Zres = Z-Zest # original minus estimated
  Zres[,colStandalone] = Z[,colStandalone] 
  
  scalest = apply(Zres[,colConnected],2,FUN=robScale,precScale=precScale)
  Zres[,colConnected] = scale(Zres[,colConnected],center=FALSE,scalest)    
  # where the generic R-function scale() scaled these residuals.
  # We don't have to scale the standalone columns, as they
  # were already standardized in the beginning.

  # Next, flag outlying cells by their large residuals:
  indcells = which(abs(Zres) > qCell) # does not flag the NAs as cells
  U[indcells] = NA
  
  } # ends the iteration
  rm(U) # we no longer use it.
  
  indcells = setdiff(indcells,col2cell(colStandalone,n=n))
  # are the indices of outlying cells in connected variables only
  indcells = unique(sort(c(indcells,indexStandalone))) 
  # are the indices of both types of outlying cells
  
  ####################################
  ##    STEP 7 : FLAGGING ROWS      ##
  ####################################      

  Ti      = integer(0)
  indrows = integer(0)
  indall  = indcells

  compT = function(rowi){ mean(pchisq(rowi^2,1),na.rm=T) - 0.5 }
   
  if(rowdetect == TRUE) {
    Ti = as.vector(apply(Zres,1,FUN=compT))    
    # calculate the test value (outlyingness) of each row:     
    Ti = scale(Ti, median(Ti,na.rm = TRUE), mad(Ti,na.rm = TRUE))
    rownames(Ti) = rownames(X)
    indrows = which(is.na(uniDetect(Ti,qCut=qRow)) & (Ti>0))
    indall  = unique(c(indcells,row2cell(indrows,n,d))) 
  } 
  
  #############################################
  ##    STEP 8: UNSTANDARDIZE AND IMPUTE     ##
  #############################################

  # compute Xest (storing in the existing matrix Zest to save space)
  Zest = sweep(sweep(Zest,2,scaleX,"*"),2,locX,"+")
  
  # compute Ximp (storing in the existing matrix X to save space)
  X[indall]  = Zest[indall]  # imputes the outlying cells of X
  X[missIndex] = Zest[missIndex] # imputes the missing values of X
  
  attr(Zres,"scaled:scale") = NULL
  attr(Ti,"scaled:center")  = NULL
  attr(Ti,"scaled:scale")   = NULL

  wnq(paste(" The computation ended at: ",date(),sep=""))
  wnq(" ")
  
  return(list(Z=Z,
              k=k,
              ngbrs=ngbrs,
              robcors=robcors,
              robslopes=robslopes,              
              Xest=Zest,
              stdResid=Zres,
              indcells=indcells,
              Ti=Ti,
              indrows=indrows,
              indall=indall,              
              Ximp=X))  
}



cellMap = function(D, R, indcells, indrows, showVals=FALSE, xlabels="", ylabels="",
                   mTitle="",xtitle="",ytitle="",xshowindex=NULL,yshowindex=NULL,
                   xblocksize=1,yblocksize=1,autolabel=TRUE,anglex=90,sizexy=1.1,
                   hjustXlabels=1,hjustYlabels=1) {
  
  # Draws a cellmap, possibly of a subset of rows and columns of the data,
  # and possibly combining cells into blocks. 
  # The inputs are:
  #
  # D            : the data matrix (required input argument)
  # R            : matrix of cell residuals (required input argument)
  # indcells     : indices of outlying cells (required input argument)
  # indrows      : indices of outlying rows (required input argument)
  # showVals     : whether to show the entries of D in the cellmap
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
  
  
  variable=rownr=rescaleoffset=NULL
  
  funcSqueeze = function(Xin,n,d,xblocksize,yblocksize) {
    # function to combine cells into blocks
    Xblock = matrix(0,nrow=n,ncol=d)
    Xblockgrad = matrix(0,nrow=n,ncol=d)
    for (i in 1:n) {
      for (j in 1:d) {
        Xsel = Xin[(1+((i-1)*yblocksize)):(i*yblocksize),
                   (1+((j-1)*xblocksize)):(j*xblocksize)]
        seltable = tabulate(Xsel,nbins =3)
        cnt0 = (yblocksize*xblocksize) - sum(seltable)
        if (sum(seltable) > 0) {
          indmax = which(seltable==max(seltable))[1]
          cntmax = seltable[indmax]
          gradmax = cntmax / (xblocksize*yblocksize)
        } else {
          indmax = 0
          gradmax = 1
        }
        Xblock[i,j] = indmax
        Xblockgrad[i,j] = gradmax
      }
    }
    return(list(X=Xblock,Xgrad=Xblockgrad))
  }
  
  n = nrow(R)
  d = ncol(R)
  
  #check input arguments 
  blockMap = FALSE
  if (xblocksize > 1 | yblocksize > 1 ){
    blockMap = TRUE
    if (xblocksize > d) stop('Input argument xblocksize cannot be larger than d')
    if (yblocksize > n) stop('Input argument yblocksize cannot be larger than n')
    if (showVals == TRUE) warning('The option showVals=TRUE cannot be combined
                                  with xblocksize or yblocksize greater than 1,
                                  so showVals is set to FALSE here.')
    showVals = FALSE
  }
  
  if(!blockMap){
    if (!all(dim(R) == dim(D))) stop('Dimensions of D and R must match')
  }
  
  if(!(blockMap & autolabel==FALSE)){
    if (length(xlabels) > 0 & length(xlabels)!= d) {
       stop(paste('Number of xlabels does not match d = ',d,sep=""))}
    if (length(ylabels) > 0 & length(ylabels)!= n) {
       stop(paste('Number of ylabels does not match n = ',n,sep=""))}
  }
  
  if(!(is.null(xshowindex) & is.null(yshowindex))){
    # here we extract the rows and columns that will be shown.
    if(is.null(xshowindex)) { xshowindex = 1:d } else {
      if(!(all(xshowindex %in% 1:d))) stop(" xshowindex goes out of bounds")}
    if(is.null(yshowindex)) { yshowindex = 1:n } else {
      if(!(all(yshowindex %in% 1:n))) stop(" yshowindex goes out of bounds")}
  
    tempMat = matrix(0,n,d)
    tempMat[indcells] = 1
    tempMat = tempMat[yshowindex,xshowindex]
    indcells = which(tempMat == 1)
  
    tempVec = rep(0,n)
    tempVec[indrows] = 1
    tempVec = tempVec[yshowindex]
    indrows = which(tempVec == 1) # also works if indrows was empty
    rm(tempMat,tempVec)
    
    if(!blockMap) D = D[yshowindex,xshowindex] 
    R = R[yshowindex,xshowindex]
    xlabels = xlabels[xshowindex]
    ylabels = ylabels[yshowindex]
    n = nrow(R)
    d = ncol(R)
  }  

  # create the matrix which indicates the color of each cell
  # 0=yellow, 1=blue, 2=red, 3=black, 4=white
  X = matrix(0,n,d)
  X[indrows,] = 3
  pcells = indcells[indcells %in% which(R>=0)]
  ncells = indcells[indcells %in% which(R<0)]
  X[ncells] = 1
  X[pcells] = 2
  if (!blockMap) X[is.na(D)] = 4
  
  if (blockMap) { # in this case the input D will be ignored  
    n = floor(n/yblocksize)
    d = floor(d/xblocksize)
    
    # create new X and Xgrad
    result = funcSqueeze(X,n,d,xblocksize,yblocksize)
    X = result$X
    Xgrad = result$Xgrad
    
    # if autolabel=F, labels{x,y} will be used for the blocks.
    if (autolabel==TRUE) { # automatically combine labels for blocks
      
      if (xblocksize>1 & length(xlabels)>0) {
        labx = xlabels
        xlabels = rep(0,d)
        for(ind in 1:d) {
          xlabels[ind] = paste(labx[(1+((ind-1)*xblocksize))],"-",
                               labx[(ind*xblocksize)],sep="")
        }
      }
      
      if (yblocksize>1 & length(ylabels)>0) {
        laby = ylabels
        ylabels = rep(0,n)
        for(ind in 1:n) {
          ylabels[ind] = paste(laby[(1+((ind-1)*yblocksize))],"-",
                               laby[(ind*yblocksize) ])
        }
      }
    } else {
      if (length(xlabels) > 0 & length(xlabels)!= d) {
        stop(paste(' autolabel=FALSE and number of xlabels is ',
                   length(xlabels),' but should be ',d,sep=""))}
      if (length(ylabels) > 0 & length(ylabels)!= n) {
        stop(paste(' autolabel=FALSE and number of ylabels is ',
                   length(ylabels),' but should be ',n,sep=""))}      
    }
    
    # Melt data matrices for cellMap
    Xdf = data.frame(cbind(seq(1,n,1),X))
    colnames(Xdf) = c("rownr",seq(1,d,1))
    rownames(Xdf) = NULL
    Xdf$rownr = with(Xdf, reorder(rownr, seq(n,1,-1)))
    mX = melt(Xdf,id.var="rownr", value.name = "CatNr") 
    
    Xgraddf = data.frame(cbind(seq(1,n,1),Xgrad))
    colnames(Xgraddf) = c("rownr",seq(1,d,1))
    rownames(Xgraddf) = NULL
    Xgraddf$rownr = with(Xgraddf, reorder(rownr, seq(n,1,-1)))
    mXgrad = melt(Xgraddf,id.var="rownr", value.name = "grad") 
    
    # Combine melted data
    mX$grad = mXgrad$grad
    mX$rescaleoffset = mXgrad$grad + 10*mX$CatNr
    scalerange = c(0,1)
    gradientends = scalerange + rep(c(0,10,20,30), each=2)
    colorends = c("yellow", "yellow", "yellow", "blue", "yellow",
                  "red", "yellow", "black")
    
  } else { # no blockMap
    
    # Melt data matrices for cellMap
    Ddf = data.frame(cbind(seq(1,n,1),D))
    colnames(Ddf) = c("rownr",colnames(D))
    rownames(Ddf) = NULL
    Ddf$rownr = with(Ddf, reorder(rownr, seq(n,1,-1)))
    mD = melt(Ddf,id.var="rownr") 
    
    Rdf = data.frame(cbind(seq(1,n,1),R))
    colnames(Rdf) = c("rownr",colnames(D))
    rownames(Rdf) = NULL
    Rdf$rownr = with(Rdf, reorder(rownr, seq(n,1,-1)))
    mR = melt(Rdf,id.var="rownr") 
    
    Xdf = data.frame(cbind(seq(1,n,1),X))
    colnames(Xdf) = c("rownr",seq(1,d,1))
    rownames(Xdf) = NULL
    Xdf$rownr = with(Xdf, reorder(rownr, seq(n,1,-1)))
    mX = melt(Xdf,id.var="rownr", value.name = "CatNr") 
    
    # Combine melted data
    mX$data = mD$value
    mX$resid = mR$value
    mX$rescaleoffset = 100*mX$CatNr
    scalerange = range(mX$resid)
    gradientends = rep(c(0,100,200,300,400), each=2)
    colorends = c("yellow", "yellow", "blue", "blue", "red", "red",
                  "black", "black", "white", "white")  
  }
  
  ylabels = rev(ylabels) # will be plotted from bottom to top
  base_size = 10
  
  ggp = ggplot(data=mX, aes(variable,rownr)) + 
    {if(blockMap) geom_tile(aes(fill=rescale(rescaleoffset,from= range(c(0,1) +
                            rep(c(0,10,20,30),each=2)) )),colour="white")} + 
    {if(!blockMap)  geom_tile(aes(fill=rescale(rescaleoffset,from= c(0,400) )),
                              colour="white")  } +
    scale_fill_gradientn(colours=colorends,values=rescale(gradientends),
                         rescaler=function(x, ...) x,oob=identity) +  
    ggtitle(mTitle) +
    coord_fixed() + 
    theme_classic(base_size = base_size*1) + 
    labs(x=xtitle,y=ytitle) + 
    scale_x_discrete(expand=c(0,0),labels=xlabels) +
    scale_y_discrete(expand=c(0,0),labels=ylabels) + 
    theme(legend.position="none",axis.ticks=element_blank(), 
          plot.title=element_text(size=base_size*2,vjust=1,face="bold"),
          axis.text.x=element_text(size=base_size*1.8,angle=anglex,
                                   hjust=hjustXlabels,vjust=0.5,colour="black"),
          axis.text.y=element_text(size=base_size*1.8,angle=0,
                                   hjust=hjustYlabels,colour="black"),
          axis.title.x=element_text(colour="black",size=base_size*sizexy,vjust=1),
          axis.title.y=element_text(colour="black",size=base_size*sizexy,vjust=0)) 
  
  if (showVals) {
    txtcol = mX$CatNr
    txtcol[txtcol==0] = "black"
    txtcol[txtcol==1] = "white"
    txtcol[txtcol==2] = "white"
    txtcol[txtcol==3] = "white"
    txtcol[txtcol==4] = "black"
    # the following line places `NA' in missing cells
    # the size of the symbols is also specified here
    ggp = ggp + geom_text(aes(label = ifelse(is.na(data),sprintf("%1.0f",data),
            round(data,1))),size = base_size*0.5, colour=txtcol, na.rm = TRUE) 
  } 
  
  return(ggp)
}


