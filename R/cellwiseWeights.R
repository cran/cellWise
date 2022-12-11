
weightedEM = function(X, w=NULL, lmin=NULL, crit=1e-4, 
                      maxiter=1000, initEst=NULL, computeloglik=F){
  #
  # Carries out a rowwise weighted EM algorithm to estimate
  # mu and Sigma of incomplete Gaussian data.
  #
  # Arguments:
  # X             : n by d data matrix or data frame.
  # w             : vector with n nonnegative rowwise (casewise)
  #                 weights. If NULL, all weights are set to 1 so
  #                 an unweighted EM is carried out.
  # lmin          : if not NULL, a lower bound on the eigenvalues 
  #                 of the estimated EM covariance matrix on the 
  #                 standardized data, to avoid singularity.  
  # crit          : convergence criterion of successive mu
  #                 and Sigma estimates.
  # maxiter       : maximal number of iteration steps.
  # initEst       : if not NULL, a list with initial estimates $mu
  #                 of the mean, $Sigma of the covariance matrix.  
  # computeloglik : if TRUE, the log(likelihood) is computed in
  #                 every step and reported. Default is FALSE
  #                 to save computation time.
  #
  # Value:
  # A list with components
  #
  # mu       : the estimated location vector.
  # Sigma    : the estimated covariance matrix.
  # impX     : the imputed data matrix.
  # niter    : the number of iteration steps taken.  
  # loglikhd : vector with the total log(likelihood) at every 
  #            iteration step. When computeloglik = FALSE this 
  #            array contains NA's.
  
  weightedStdev = function(x,w=rep(1,length(x))){
    # Auxiliary function.
    n = length(x)
    if(length(w) != n)    stop(" length(weights) differs from n")
    w = w[!is.na(x)]
    x = x[!is.na(x)]
    if(sum(is.na(w)) > 0) stop(" there are NA weights")
    if(sum(w < 0) > 0)    stop(" there are negative weights")
    sumw  = sum(w)
    if(sumw == 0)         stop(" all weights are zero")
    wmean = sum(w*x)/sumw
    wvar  = sum(w*(x-wmean)*(x-wmean))/sumw
    sqrt(wvar)
  }
  
  X = as.matrix(X)
  if(sum(as.vector(is.infinite(X)))>0)
    stop(" X has infinite values")
  # then do e.g.:  X[is.infinite(X)] = NA
  n = dim(X)[1]
  d = dim(X)[2] 
  if(d < 2) stop(" X must have at least 2 columns")
  cnames = colnames(X)
  # First check the weights:
  if(is.null(w)) { w = rep(1,n) } else {
    w = as.vector(w)
    if(length(w) != n)    stop(" length(weights) differs from n")  
    if(sum(is.na(w)) > 0) stop(" there are NA weights")  
    if(sum(w < 0) > 0)    stop(" there are negative weights") 
    if(sum(w==0) > 0)     warning(paste(c(" there are case weights",
                                          " of 0, so the imputed data will have fewer rows"),sep=""))
    X = X[w>0,] # only keep rows with w>0. Preserves row names.
    w = w[w>0]  # only keep rows with w>0.
  }
  # Now check the (possibly reduced) dataset X:
  if(any(colSums(!is.na(X))==0))
    stop(" X has columns with only NA's")
  nmiss = rowSums(!is.na(X)) 
  if(any(nmiss==0)) {
    warning(paste(c(" X had rows with only NA's, so the imputed",
                    " data will have fewer rows"),sep=""))
    X = X[nmiss!=0,] # only keep non-missing rows. Preserves row names.
    w = w[nmiss!=0]    
  }
  n = dim(X)[1]   # may be smaller than the original n  
  if(n <= d) stop(" n should be larger than d")
  wnorm = sum(w)/n
  w = w/wnorm     # a stable range for w, for numerical accuracy.
  sumw = sum(w)   # compute this once
  sqrtw = sqrt(w) # compute this once
  #
  # We now center and scale the columns of X, and
  # will transform back at the end.
  #
  mu0 = apply(X,2,FUN=weighted.mean,w=w,na.rm=T)
  sc0 = apply(X,2,FUN=weightedStdev,w=w)    
  if(sum(sc0==0) > 0) stop(paste(c(" X has at least one variable",
                                   " with zero scale",sep="")))   
  diagSc0 = diag(sc0)
  z       = t((t(X)-mu0)/sc0)
  #
  # Extract missingness patterns:
  Rmispat = !unique(is.na(z)) # all different _R_ow missingness patterns
  mode(Rmispat) = 'numeric'   # NA=1, observed=0
  #
  patZ = !is.na(z)
  mode(patZ) = 'numeric' # here NA=0, observed=1
  patR = apply(Rmispat,1,paste,collapse='')
  # turns row patterns into strings [example: c(1,0,1,1) -> '1011']
  patZ = apply(patZ,1,paste,collapse='')
  Imispat = list()
  Omispat = list()
  Mmispat = list()
  Nmispat = dim(Rmispat)[1] # _N_umber of different missingness patterns
  for (s in 1:Nmispat){
    Imispat[[s]] = which(patZ==patR[s])  # row _I_ndices with pattern s
    Omispat[[s]] = which(Rmispat[s,]==1) # _O_bserved cols in pattern s
    Mmispat[[s]] = which(Rmispat[s,]==0) # _M_issing cols in pattern s
  }
  #
  # Initial estimates:
  #
  impZ = z           # will hold the imputed matrix
  impZ[is.na(z)] = 0 # NA's are first imputed by zeroes
  # Zeroes make sense because we centered the data above.
  # Next, compute the observed part of the statistic (T1,T2)
  # where T1 is the sum of the rows and T2 is the sscp matrix:
  T1obs = colSums(w*impZ)
  T2obs = t(sqrtw*impZ)%*%(sqrtw*impZ)
  if(is.null(initEst)){
    mu       = T1obs/sumw  
    Sigma    = ((1/sumw)*T2obs)-mu%*%t(mu)
    invSigma = mpinv(Sigma)$Inv  
  } else { # added this part
    mu       = (initEst$mu-mu0)/sc0
    Sigma    = diag(1/sc0) %*% initEst$Sigma %*% diag(1/sc0)
    Sigma    = truncEig(Sigma, lmin=1e-6)
    invSigma = diagSc0 %*% mpinv(Sigma)$Inv %*% diagSc0
  }
  #
  # Prepare for iterations:
  #
  param    = c(mu,as.vector(invSigma))
  imperror = rep(0,maxiter)
  err      = Inf
  iter     = 0
  allMu    = matrix(0,nrow=d,ncol=maxiter)
  allSigma = array(0,dim=c(d,d,maxiter))
  loglikhd = rep(NA,maxiter)
  #
  ## Iterations:
  #
  while ((err>crit)&(iter<maxiter)){
    param.old = param
    impZ.old  = impZ
    T1 = T1obs
    T2 = T2obs
    for (s in 1:Nmispat){
      mis = Mmispat[[s]]   # extract M(s) # length(mis) <= d
      obs = Omispat[[s]]   # extract O(s) # length(obs) <= d
      imisp = Imispat[[s]] # extract I(s) # length(imisp) <= n
      wmispat = w[imisp]   # extract weights of these rows
      if (length(mis)>0){ # i.e. skip the fully observed rows
        #
        ## Part 1: compute parameters of z(mis)|z(obs)
        # kmisinv = solve(invSigma[mis,mis]) # replaced by mpinv:
        kmisinv = mpinv(invSigma[mis,mis])$Inv
        B       = invSigma[obs,mis,drop=FALSE]%*%kmisinv # regr coeffs
        #
        ## Part 2: E-step: replace the missings by their expectation
        mumis = matrix(mu[mis],length(imisp),length(mis),byrow=T)
        muobs = matrix(mu[obs],length(imisp),length(obs),byrow=T)
        zmishat = mumis - (z[imisp,obs] - muobs)%*%B
        impZ[imisp,mis] = zmishat # imputed data
        #
        ## Part 3: update the statistics T1 and T2
        T1[mis]     = T1[mis]+colSums(wmispat*zmishat)
        T2[obs,mis] = T2[obs,mis]+
          t(wmispat*z[imisp,obs,drop=FALSE])%*%zmishat
        T2[mis,obs] = t(T2[obs,mis])
        tempmat     = sqrt(wmispat)*zmishat
        T2[mis,mis] = T2[mis,mis] +
          t(tempmat)%*%tempmat + sum(wmispat)*kmisinv
      }
    }
    ## M-step: apply maximum likelihood
    mu    = T1/sumw                  # weighted MLE location
    Sigma = ((1/sumw)*T2)-mu%*%t(mu) # weighted MLE scatter matrix
    if(!is.null(lmin)){
      Sigma  = truncEig(Sigma, lmin = lmin, lmax = NULL)
    }
    invSigma = mpinv(Sigma)$Inv
    param = c(mu,as.vector(invSigma))
    err = sqrt(sum((param-param.old)^{2}))/(1+sqrt(sum(param^{2})))
    iter = iter+1
    imperror[iter]   = sum((impZ-impZ.old)^{2})/sum(impZ^{2})
    tempMu           = mu0 + mu*sc0
    allMu[,iter]     = tempMu
    tempSigma        = diagSc0 %*% Sigma %*% diagSc0
    allSigma[,,iter] = tempSigma
    loglik = 0
    if(computeloglik == T){
      for (s in 1:Nmispat){
        tempObs  = Omispat[[s]] # set of observed columns
        subMu    = tempMu[tempObs]
        subSigma = tempSigma[tempObs,tempObs,drop=F]
        invSubSigma = mpinv(subSigma)$Inv 
        imisp = Imispat[[s]] # extract the rows with pattern s
        ns = sum(w[imisp])   # total weight of this pattern
        ds = length(tempObs) # number of observed columns
        Ms = 0
        for (i in imisp){
          Ms = Ms + w[i]*(X[i,tempObs]-subMu)%*%t(X[i,tempObs]-subMu)
        }
        loglik = loglik-((ns*ds/2)*log(2*pi))-
          ((ns/2)*log(det(subSigma)))-((sum(diag(invSubSigma%*%Ms)))/2)
      }
      loglikhd[iter] = wnorm*loglik # scales back to original weights
    } # ends computation of loglikelihood
  } # ends while loop of iterations
  if (err>crit) warning(
    " reached maximal number of iterations without attaining crit")  
  mu = allMu[,iter]
  names(mu) = cnames
  Sigma = allSigma[,,iter]
  rownames(Sigma) = colnames(Sigma) = cnames
  # Sigmai = diag(1/sc0) %*% invSigma %*% diag(1/sc0)
  impX  = X
  impX2 = t(t(impZ)*sc0 + mu0)
  impX[is.na(X)] = impX2[is.na(X)]
  loglikhd = loglikhd[1:iter]
  list(mu=mu,Sigma=Sigma,impX=impX,niter=iter,loglikhd=loglikhd)
}  


unpack = function(X,W){
  # Unpacks cellwise weighted data.
  #
  # This function transforms a dataset X with cellwise weights W
  # to an extended data matrix U with more rows, and containing 
  # more NA's. Its rows have the case weights v.
  #
  # Arguments:
  # X : An n by d data matrix or data frame. Must be given.
  #     X is allowed to contain NA's.
  # W : An n by d matrix of nonnegative cellwise weights.
  #     Must be given. W is not allowed to contain NA's.
  #
  # Value:
  # A list with components
  #
  # U : unpacked data matrix, with the same columns
  #     as X but typically more rows.
  # v : vector with the rowwise (=casewise) weights of U.
  
  # First check the input:
  X = as.matrix(X)
  n = nrow(X)
  d = ncol(X)
  if(sum(as.vector(is.infinite(X[,-1])))>0)
    stop(" X has infinite values")
  # then do e.g.:  X[is.infinite(X)] = NA  
  if(d < 2) stop(" X must have at least 2 columns")
  W = as.matrix(W)
  # First check the weights:
  if(nrow(W) != n) stop(" nrow(W) differs from nrow(X)")  
  if(ncol(W) != d) stop(" ncol(W) differs from ncol(X)")
  if(sum(as.vector(is.infinite(W)))>0)
    stop(" W has infinite values")  
  W[is.na(X)] = 0 # NA's get zero weight
  X[is.na(X)] = 0 # from here on X has no NA's  
  if(sum(is.na(W)) > 0) stop(" there are NA weights")  
  if(sum(W < 0) > 0)    stop(" there are negative weights")
  if(is.null(rownames(X))) rownames(X) = seq_len(n)
  rnames = rownames(X) # store row names, to put back later
  if(is.null(colnames(X))) colnames(X) = seq_len(d)
  cnames = colnames(X) # store column names too
  X = cbind(seq_len(n),X) # new column 1 is original numbering of rows
  wnz = rowSums(W>0)  
  if(any(wnz == 0)){
    X = X[wnz>0,] # here rows of X can be dropped
    W = W[wnz>0,]        
    warning(paste(c(" There were rows with only zero weights,",
                    " we dropped them from both X and W"),sep="")) 
  }
  wnz = colSums(W>0)  
  if(any(wnz==0)) stop(" There are columns with only zero weights")  
  colRanges = rep(0,d)
  for(j in 1:d){
    Xj = X[,j+1] # because we added a "fake" first column to X.
    Wj = W[,j]
    Xj = Xj[Wj > 0]
    if(length(Xj)>0) colRanges[j] = max(Xj)-min(Xj)
  }
  if(any(colRanges==0))
    stop(" X has columns with constant relevant entries")  
  #
  # Transform (X,W) to (U,v):
  #
  XW = cbind(X, W) # column 1 is the original row number.
  #
  unpackRow = function(XWi){ # internal function
    # XWi will be row i of XW = cbind(X,W).
    num = XWi[1]  # the row number in the original X matrix.
    XWi = XWi[-1] # take it out again.
    d = length(XWi)/2 # We checked d before applying unpackRow(). 
    Xi  = XWi[1:d]         # a row of the original X matrix.
    Wi  = XWi[(d+1):(2*d)] # corresponding row of w matrix.
    out = sort(Wi,decreasing=T,index.return=T)
    uniw = unique(out$x) # the unique weights in this row Wi.
    if(uniw[length(uniw)]>0) uniw=c(uniw,0)
    ontbr = d - length(uniw) + 1 # we want d+1 entries
    if(ontbr >= 0) uniw = c(uniw, rep(uniw[length(uniw)], ontbr))
    # This adds fake values that will yield zeroes in v.
    NArow = rep(NA,d) # initialize as NA 
    vu = integer(0)   # initialize as empty
    for(j in seq_len(d)){ # all outputs vu will have same length.
      uj = NArow
      indx = (Wi >= uniw[j])
      uj[indx] = Xi[indx]
      vj = uniw[j] - uniw[j+1]
      vu = c(vu, num, vj, uj) # each time we append d+2 entries,
      # so the output of unpackRow is a vector with length d*(d+2).
    }
    vu
  }
  #
  U = apply(XW, 1, unpackRow) # has n*d*(d+2) entries
  U = matrix(as.vector(U), ncol=d+2, byrow=T)
  # This gives U the right dimensions.
  if(!is.null(rnames)) rownames(U) = rnames[U[,1]] 
  # puts back the correct row names of X. 
  U = U[,-1] # removes the "fake" first column.
  U = U[U[,1] > 0,] # remove rows with weight 0
  v = U[,1]  # splits off the weights v
  U = U[,-1] # the unpacked matrix U
  colnames(U) = cnames
  list(U=U,v=v)  
}


cwLocScat = function(X, W, methods = "all", lmin = 1e-3,
                     crit = 1e-12, maxiter= 1000, 
                     initCwCov = FALSE, initEst = NULL){
  #
  # Computes different estimators of multivariate location
  # and scatter for cellwise weighted data.
  #
  # Arguments:
  # X         : An n by d data matrix or data frame. Must be
  #             given. X is allowed to contain NA's.
  # W         : An n by d matrix of nonnegative cellwise weights.
  #             Must be given. W is not allowed to contain NA's.  
  # methods   : either "all" or "explicit". If "explicit"
  #             only the explicit estimates cwMean, cwCov
  #             and sqrtCov are computed. If "all" (the
  #             default) also the cellwise MLE is carried out,
  #             yielding cwMLEmu and cwMLEsigma.  
  # lmin      : if not NULL, a lower bound on the eigenvalues 
  #             of the estimated covariance matrices on the 
  #             standardized data, to avoid singularity.
  # crit      : convergence criterion of successive mu
  #             and Sigma estimates in the EM algorithm.
  # maxiter   : maximal number of iteration steps in EM.
  # initcwCov : if TRUE, uses the weighted mean and cwCov
  #             as initial estimates for the weighted EM.
  # initEst   : if not NULL, a list with initial estimates $mu
  #             of the mean, $Sigma of the covariance matrix,
  #             for the weighted EM. Has no effect when
  #             initCwCov = TRUE.
  #   
  # Outputs:
  # cwMean     : the explicit cellwise weighted mean.
  # cwCov      : explicit cellwise weighted covariance matrix.
  #              Is asymptotically normal but not necessarily
  #              PSD (unless a nonnegative lmin was specified).
  # sqrtCov    : the cellwise weighted covariance matrix of Van 
  #              Aelst et al (2011). Also asymptotically normal 
  #              but not necessarily PSD (unless a nonnegative 
  #              lmin was specified).
  # cwMLEmu    : the location estimate obtained by the cwMLE.
  # cwMLEsigma : the covariance matrix obtained by the cwMLE.
  #              Is PSD when the EM algorithm converges.
  
  weightedStdev = function(x,w=rep(1,length(x))){
    # Auxiliary function.
    n = length(x)
    if(length(w) != n)    stop(" length(weights) differs from n")
    w = w[!is.na(x)]
    x = x[!is.na(x)]
    if(sum(is.na(w)) > 0) stop(" there are NA weights")
    if(sum(w < 0) > 0)    stop(" there are negative weights")
    sumw  = sum(w)
    if(sumw == 0)         stop(" all weights are zero")
    wmean = sum(w*x)/sumw
    wvar  = sum(w*(x-wmean)*(x-wmean))/sumw
    sqrt(wvar)
  }
  
  if(!(methods %in% c("explicit","all"))){
    stop('methods should be either "explicit" or "all".') }
  X = as.matrix(X)
  d = ncol(X)
  W = as.matrix(W)
  W[is.na(X)] = 0 # NA's get zero weight
  #
  # The next call will also check the inputs:
  out = unpack(X,W)
  U = out$U
  v = out$v
  cnames = colnames(U)
  #
  # Initial standardization:
  #
  mu0 = apply(U,2,FUN=weighted.mean,w=v,na.rm=T)
  # This is the cellwise weighted mean.
  sc0 = apply(U,2,FUN=weightedStdev,w=v)  
  if(sum(sc0==0) > 0) stop(paste(c(" X has at least one variable",
                                   " with zero scale",sep="")))   
  diagSc0 = diag(sc0)
  X       = t((t(X)-mu0)/sc0) # reuse name X to save space
  #
  # Compute cwCov: explicit cellwise weighted covariance matrix
  #
  cwCov = matrix(rep(0,d*d),nrow=d)
  for(j in 1:d){
    for(k in 1:d){
      minweights = pmin(W[,j],W[,k])
      covterms = minweights*X[,j]*X[,k]
      indobs = which(!is.na(covterms))
      numer = sum(covterms[indobs])
      denom = sum(minweights[indobs])
      if(denom <= 0) { cwCov[j,k] = 0 } else {
        cwCov[j,k] = numer/denom
      }
    }
  }
  cwCov = truncEig(cwCov,lmin=lmin,lmax = NULL)
  cwCov = diagSc0 %*% cwCov %*% diagSc0
  rownames(cwCov) = colnames(cwCov) = cnames
  #
  # Compute sqrtCov: from Van Aelst et al (2011).
  #
  sqrtW = sqrt(W)
  sqrtCov = matrix(rep(0,d*d),nrow=d)
  for(j in 1:d){
    for(k in 1:d){
      prodweights = sqrtW[,j]*sqrtW[,k]
      covterms = prodweights*X[,j]*X[,k]
      indobs = which(!is.na(covterms))
      numer = sum(covterms[indobs])
      denom = sum(prodweights[indobs])
      if(denom <= 0) { sqrtCov[j,k] = 0 } else {
        sqrtCov[j,k] = numer/denom
      }
    }
  }
  sqrtCov = truncEig(sqrtCov,lmin=lmin,lmax=NULL)
  sqrtCov = diagSc0 %*% sqrtCov %*% diagSc0
  rownames(sqrtCov) = colnames(sqrtCov) = cnames
  #
  # Compute weighted MLE estimates:
  #
  if(methods == "all"){  
    if(initCwCov == TRUE) initEst = list(mu=mu0, Sigma=cwCov)
    EMfit = weightedEM(U, v, crit=crit, maxiter=maxiter,
                       computeloglik=F, initEst=initEst, lmin=lmin)
    res = list(cwMean = mu0, cwCov = cwCov, sqrtCov = sqrtCov,
               cwMLEmu = EMfit$mu, cwMLEsigma = EMfit$Sigma,
               cwMLEiter = EMfit$niter)
  } else {
    res = list(cwMean = mu0, cwCov = cwCov, sqrtCov = sqrtCov)
  }
  return(res)
}
