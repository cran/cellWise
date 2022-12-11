



.detmcd_fixedcenter <- function(x, h, hsets.init = NULL,
                                save.hsets = missing(hsets.init), 
                                full.h = save.hsets, scalefn, 
                                maxcsteps = 200,
                                warn.nonconv.csteps = getOption("robustbase:warn.nonconv.csteps",TRUE),
                                warn.wrong.obj.conv = getOption("robustbase:warn.wrong.obj.conv",FALSE),
                                trace = as.integer(trace), names = TRUE) {
  # function based on robustbase:::.detmcd
  # but now with the option fixedCenter added:
  fixedCenter <- TRUE # JR
  
  stopifnot(length(dx <- dim(x)) == 2, h == as.integer(h), 
            h >= 1)
  n <- dx[1]
  p <- dx[2]
  stopifnot(p >= 1, n >= 1)
  
  if (fixedCenter) {# JR
    scalefn <- function(x) {1.4826 * median(abs(x))} # MAD around zero
    ## Maybe we should also allow other scale functions?
    vnms <- colnames(x)
    z <- doScale(unname(x), center = rep(0, p), scale = scalefn)
    z.center <- z$center
    z.scale <- z$scale
    z <- z$x # JR
  } else {
    scalefn <- robScalefn(scalefn, n)
    vnms <- colnames(x)
    z <- doScale(unname(x), center = median, scale = scalefn)
    z.center <- z$center
    z.scale <- z$scale
    z <- z$x
  }
  
  if (is.null(hsets.init)) { # then use DetMCD starts
    # for scaled = TRUE, r6pack does no additional centering/scaling,
    # so this is fine for both fixedCenter TRUE and FALSE: # JR
    hsets.init <- robustbase::r6pack(z, h = h, full.h = full.h,
                                     scaled = TRUE, 
                                     scalefn = scalefn)
    dh <- dim(hsets.init)
  }
  else { # then use the _given_ h'-subsets
    if (is.vector(hsets.init)) 
      hsets.init <- as.matrix(hsets.init)
    dh <- dim(hsets.init)
    if (dh[1] < h || dh[2] < 1) 
      stop("'hsets.init' must be a  h' x L  matrix (h' >= h) of observation indices")
    if (full.h && dh[1] != n) 
      warning("'full.h' is true, but 'hsets.init' has less than n rows")
    if (min(hsets.init) < 1 || max(hsets.init) > n) 
      stop("'hsets.init' must be in {1,2,...,n}; n = ", 
           n)
  }
  nsets <- ncol(hsets.init) # = number of h'-sets
  hset.csteps <- integer(nsets) # = rep(0,nsets) # initialization
  bestobj <- Inf
  for (i in 1:nsets) { # loop over the initial h'-sets
    if (trace) {
      if (trace >= 2) 
        cat(sprintf("H-subset %d = observations c(%s):\n-----------\n", 
                    i, pasteK(hsets.init[1:h, i])))
      else cat(sprintf("H-subset %d: ", i))
    }
    for (j in 1:maxcsteps) { # iterate C-steps
      if (j == 1) {
        obs_in_set <- hsets.init[1:h, i]
      }
      else {
        score <- (z - rep(svd$center, each = n)) %*% svd$loadings
        mah <- robustbase::mahalanobisD(score, center = FALSE, sd = sqrt(abs(svd$eigenvalues)))
        obs_in_set <- sort.list(mah)[1:h]
      }
      # C-step is done by SVD in robustbase:::classPC :
      svd <- robustbase::classPC(z[obs_in_set, , drop = FALSE],
                                 center = !fixedCenter,
                                 signflip = FALSE,
                                 via.svd = TRUE) # JR
      if (fixedCenter) { svd$center <- rep(0, p) }# JR
      obj <- sum(log(svd$eigenvalues)) # = log(det)
      if (svd$rank < p) {
        stop("More than h of the observations lie on a hyperplane.")
        exactfit <- TRUE
      }
      if (j >= 2 && obj == prevdet) {
        if (identical(obs_in_set, prevobs)) 
          break
        if (warn.wrong.obj.conv) 
          warning(sprintf("original detmcd() wrongly declared c-step convergence (obj=%g, i=%d, j=%d)", 
                          obj, i, j))
      }
      prevdet <- obj
      prevobs <- obs_in_set
    }
    hset.csteps[i] <- j
    if (trace) 
      cat(sprintf("%3d csteps, obj=log(det|.|)=%g", 
                  j, obj))
    if (obj < bestobj) {
      if (trace) 
        cat(" = new optim.\n")
      bestset <- obs_in_set
      bestobj <- obj
      initmean <- svd$center
      L <- svd$loadings
      initcov <- tcrossprod(L * rep(svd$eigenvalues, each = nrow(L)), 
                            L)
      ind.best <- i
    }
    else if (obj == bestobj) 
      ind.best <- c(ind.best, i)
    else if (trace) 
      cat("\n")
  }
  if (warn.nonconv.csteps && any(eq <- hset.csteps == maxcsteps)) {
    p1 <- paste(ngettext(sum(eq), "Initial set", "Initial sets"), 
                pasteK(which(eq)))
    warning(sprintf("%s did not converge in maxcsteps=%d concentration steps", 
                    p1, maxcsteps), domain = NA)
  }
  raw.cov <- initcov * tcrossprod(z.scale)
  raw.center <- initmean * z.scale + z.center
  if (names) {
    dimnames(raw.cov) <- list(vnms, vnms)
    names(raw.center) <- vnms
  }
  raw.objective <- bestobj + 2 * sum(log(z.scale))
  list(initcovariance = raw.cov, initmean = raw.center, best = bestset, 
       mcdestimate = raw.objective, iBest = ind.best, n.csteps = hset.csteps, 
       initHsets = if (save.hsets) hsets.init, exactfit = 0)
}


covMcd2 <- function (x, cor = FALSE, raw.only = FALSE, center = NULL, 
                     alpha = control$alpha,
                     nsamp = control$nsamp, nmini = control$nmini, 
                     kmini = control$kmini, scalefn = control$scalefn,
                     maxcsteps = control$maxcsteps, initHsets = NULL,
                     save.hsets = FALSE, names = TRUE, seed = control$seed,
                     tolSolve = control$tolSolve, trace = control$trace, 
                     use.correction = control$use.correction,
                     wgtFUN = control$wgtFUN, 
                     control = rrcov.control()) 
{
  # code based on robustbase::covMCD.
  # added option fixedCenter, only for detMCD:
  if(is.null(center)){ fixedCenter <- FALSE
  } else {
    fixedCenter <- TRUE # JR
    givenCenter <- as.vector(center) # PR
    # if(!is.vector(center)) stop("center must be a vector") # PR
    nsamp <-  "deterministic" # JR
  }
  logdet.Lrg <- 50
  if (length(seed) > 0) {
    if (length(seed) < 3 || seed[1L] < 100) 
      stop("invalid 'seed'. Must be compatible with .Random.seed !")
    if (exists(".Random.seed", envir = .GlobalEnv, 
               inherits = FALSE)) {
      seed.keep <- get(".Random.seed", envir = .GlobalEnv, 
                       inherits = FALSE)
      on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
    }
    assign(".Random.seed", seed, envir = .GlobalEnv)
  }
  defCtrl <- if (missing(control)) 
    control
  else rrcov.control()
  if (missing(wgtFUN)) 
    getDefCtrl("wgtFUN", defCtrl)
  if (is.null(nmini)) 
    getDefCtrl("nmini", defCtrl)
  if (is.numeric(nsamp) && nsamp <= 0) 
    stop("Invalid number of trials nsamp = ", nsamp, 
         "!")
  if (is.data.frame(x)) 
    x <- data.matrix(x, rownames.force = FALSE)
  else if (!is.matrix(x)) 
    x <- matrix(x, length(x), 1, dimnames = if (names) 
      list(names(x), deparse(substitute(x))))
  if (!names) 
    dimnames(x) <- NULL
  ok <- is.finite(x %*% rep.int(1, ncol(x)))
  x <- x[ok, , drop = FALSE]
  if (!length(dx <- dim(x))) 
    stop("All observations have missing values!")
  n <- dx[1]
  p <- dx[2]
  if (names) 
    dimn <- dimnames(x)
  if(fixedCenter){ # PR
    pc = length(givenCenter)
    if(pc != p) stop(paste0(
      "length(center) = ",pc," should equal ncol(x) = ",p))
    savedx <- x # to put in final output
    x <- scale(x, center = givenCenter, scale = FALSE)
  }
  h <- robustbase::h.alpha.n(alpha, n, p)
  if (n <= p + 1) 
    stop(if (n <= p) 
      "n <= p -- you can't be serious!"
      else "n == p+1  is too small sample size for MCD")
  if (n < 2 * p) {
    warning("n < 2 * p, i.e., possibly too small sample size")
  }
  if (h > n) 
    stop("Sample size n  <  h(alpha; n,p) := size of \"good\" subsample")
  else if (2 * h < n) 
    warning("subsample size\t h < n/2  may be too small")
  if (is.character(wgtFUN)) {
    if (is.function(mkWfun <- .wgtFUN.covMcd[[wgtFUN]])) 
      wgtFUN <- mkWfun(p = p, n = n, control)
  }
  if (!is.function(wgtFUN)) 
    stop(gettextf("'wgtFUN' must be a function or one of the strings %s.", 
                  pasteK(paste0("\"", names(.wgtFUN.covMcd), 
                                "\""))), domain = NA)
  raw.cnp2 <- cnp2 <- c(1, 1)
  ans <- list(call = match.call(), nsamp = nsamp, method = sprintf("MCD(alpha=%g ==> h=%d)", 
                                                                   alpha, h))
  if (h == n) {
    if (fixedCenter) {# JR
      mcd <- (t(x) %*% x)/n
      loc <- rep(0, p)
    } else {
      mcd <- cov(x)
      loc <- as.vector(colMeans(x))
    }
    obj <- determinant(mcd, logarithm = TRUE)$modulus[1]
    if (-obj/p > logdet.Lrg) {
      ans$cov <- mcd
      if (names) 
        dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
      if (cor) 
        ans$cor <- cov2cor(ans$cov)
      ans$center <- loc
      if (names && length(dimn[[2]])) 
        names(ans$center) <- dimn[[2]]
      ans$n.obs <- n
      ans$singularity <- list(kind = "classical")
      weights <- 1
    }
    else {
      mah <- mahalanobis(x, loc, mcd, tol = tolSolve)
      weights <- wgtFUN(mah)
      sum.w <- sum(weights)
      ans <- c(ans, cov.wt(x, wt = weights, cor = cor, 
                           center = !fixedCenter))# JR
      if (fixedCenter) {ans$center <- rep(0, p)}# JR
      if (sum.w != n) {
        cnp2[1] <- .MCDcons(p, sum.w/n)
        ans$cov <- ans$cov * cnp2[1]
      }
      obj <- determinant(mcd, logarithm = TRUE)$modulus[1]
      if (-obj/p > logdet.Lrg) {
        ans$singularity <- list(kind = "reweighted.MCD")
      }
      else {
        mah <- mahalanobis(x, ans$center, ans$cov, tol = tolSolve)
        weights <- wgtFUN(mah)
      }
    }
    ans$alpha <- alpha
    ans$quan <- h
    ans$raw.cov <- mcd
    ans$raw.center <- loc
    if (names && !is.null(nms <- dimn[[2]])) {
      names(ans$raw.center) <- nms
      dimnames(ans$raw.cov) <- list(nms, nms)
    }
    ans$crit <- obj
    ans$method <- paste(ans$method, "\nalpha = 1: The minimum covariance determinant estimates based on", 
                        n, "observations \nare equal to the classical estimates.")
    ans$mcd.wt <- rep.int(NA, length(ok))
    ans$mcd.wt[ok] <- weights
    if (names && length(dimn[[1]])) 
      names(ans$mcd.wt) <- dimn[[1]]
    ans$wt <- NULL
    ans$X <- x
    if (names) {
      if (length(dimn[[1]])) 
        dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
      else dimnames(ans$X) <- list(seq(along = ok)[ok], 
                                   NULL)
    }
    if (trace) 
      cat(ans$method, "\n")
    ans$raw.cnp2 <- raw.cnp2
    ans$cnp2 <- cnp2
    class(ans) <- "mcd"
    return(ans)
  }
  mcd <- if (nsamp == "deterministic") {
    ans$method <- paste("Deterministic", ans$method)
    if (fixedCenter) {# JR
      .detmcd_fixedcenter(x, h, hsets.init = initHsets, save.hsets = save.hsets, 
                          scalefn = scalefn, maxcsteps = maxcsteps, trace = as.integer(trace), 
                          names = names)# JR
      
    } else {
      .detmcd(x, h, hsets.init = initHsets, save.hsets = save.hsets, 
              scalefn = scalefn, maxcsteps = maxcsteps, trace = as.integer(trace), 
              names = names)
    }
  }
  else {
    ans$method <- paste0("Fast ", ans$method, "; nsamp = ", 
                         nsamp, "; (n,k)mini = (", nmini, ",", 
                         kmini, ")")
    .fastmcd(x, h, nsamp, nmini, kmini, trace = as.integer(trace))
  }
  calpha <- .MCDcons(p, h/n)
  correct <- if (use.correction) 
    .MCDcnp2(p, n, alpha)
  else 1
  raw.cnp2 <- c(calpha, correct)
  if (p == 1) {
    ans$method <- paste("Univariate", ans$method)
    scale <- sqrt(calpha * correct) * as.double(mcd$initcovariance)
    center <- as.double(mcd$initmean)
    if (abs(scale - 0) < 1e-07) {
      ans$singularity <- list(kind = "identicalObs", 
                              q = h)
      ans$raw.cov <- ans$cov <- matrix(0)
      ans$raw.center <- ans$center <- center
      ans$n.obs <- n
      ans$alpha <- alpha
      ans$quan <- h
      if (names && !is.null(nms <- dimn[[2]][1])) {
        names(ans$raw.center) <- names(ans$center) <- nms
        dimnames(ans$raw.cov) <- dimnames(ans$cov) <- list(nms, 
                                                           nms)
      }
      ans$crit <- -Inf
      weights <- as.numeric(abs(x - center) < 1e-07)
    }
    else {
      weights <- wgtFUN(((x - center)/scale)^2)
      sum.w <- sum(weights)
      ans <- c(ans, cov.wt(x, wt = weights, cor = cor, 
                           center = !fixedCenter)) # JR
      if (fixedCenter) {ans$center <- rep(0, p)} # JR
      if (sum.w != n) {
        cdelta.rew <- .MCDcons(p, sum.w/n)
        correct.rew <- if (use.correction) 
          .MCDcnp2.rew(p, n, alpha)
        else 1
        cnp2 <- c(cdelta.rew, correct.rew)
        ans$cov <- cdelta.rew * correct.rew * ans$cov
      }
      ans$alpha <- alpha
      ans$quan <- h
      ans$raw.cov <- as.matrix(scale^2)
      ans$raw.center <- as.vector(center)
      if (names && !is.null(nms <- dimn[[2]][1])) {
        dimnames(ans$raw.cov) <- list(nms, nms)
        names(ans$raw.center) <- nms
      }
      ans$crit <- log(sum(sort((x - as.double(mcd$initmean))^2, 
                               partial = h)[1:h])/max(1, h - 1))
      center <- ans$center
      scale <- as.vector(sqrt(ans$cov))
      weights <- wgtFUN(((x - center)/scale)^2)
    }
  }
  else {
    mcd$initcovariance <- matrix(calpha * correct * mcd$initcovariance, 
                                 p, p)
    if (raw.only || mcd$exactfit != 0) {
      dim(mcd$coeff) <- c(5, p)
      ans$cov <- ans$raw.cov <- mcd$initcovariance
      ans$center <- ans$raw.center <- as.vector(mcd$initmean)
      if (names && !is.null(nms <- dimn[[2]])) {
        dimnames(ans$cov) <- list(nms, nms)
        names(ans$center) <- nms
      }
      ans$n.obs <- n
      if (raw.only) {
        ans$raw.only <- TRUE
      }
      else {
        if (!(mcd$exactfit %in% c(1, 2, 3))) 
          stop("Unexpected 'exactfit' code ", mcd$exactfit, 
               ". Please report!")
        ans$singularity <- list(kind = "on.hyperplane", 
                                exactCode = mcd$exactfit, p = p, h = h, count = mcd$kount, 
                                coeff = mcd$coeff[1, ])
      }
      ans$alpha <- alpha
      ans$quan <- h
      if (names && !is.null(nms <- dimn[[2]])) {
        names(ans$raw.center) <- nms
        dimnames(ans$raw.cov) <- list(nms, nms)
      }
      ans$crit <- -Inf
      weights <- mcd$weights
    }
    else {
      mah <- mahalanobis(x, mcd$initmean, mcd$initcovariance, 
                         tol = tolSolve)
      weights <- wgtFUN(mah)
      sum.w <- sum(weights)
      ans <- c(ans, cov.wt(x, wt = weights, cor = cor, 
                           center = !fixedCenter))# JR
      if (fixedCenter) {ans$center <- rep(0, p)}# JR
      
      sing.rewt <- any(apply(ans$cov == 0, 2, all))
      if (!sing.rewt && sum.w != n) {
        cdelta.rew <- .MCDcons(p, sum.w/n)
        correct.rew <- if (use.correction) 
          .MCDcnp2.rew(p, n, alpha)
        else 1
        cnp2 <- c(cdelta.rew, correct.rew)
        ans$cov <- cdelta.rew * correct.rew * ans$cov
      }
      ans$best  <- sort(as.vector(mcd$best))
      ans$alpha <- alpha
      ans$quan  <- h
      ans$raw.cov <- mcd$initcovariance
      ans$raw.center <- as.vector(mcd$initmean)
      if (names && !is.null(nms <- dimn[[2]])) {
        names(ans$raw.center) <- nms
        dimnames(ans$raw.cov) <- list(nms, nms)
      }
      ans$raw.weights <- weights
      ans$crit        <- mcd$mcdestimate
      ans$raw.mah     <- mah
      if (sing.rewt || -determinant(ans$cov, logarithm = TRUE)$modulus[1]/p > 
          logdet.Lrg) {
        ans$singularity <- list(kind = paste0("reweighted.MCD", 
                                              if (sing.rewt) "(zero col.)"))
        ans$mah <- mah
      }
      else {
        mah <- mahalanobis(x, ans$center, ans$cov, tol = tolSolve)
        ans$mah <- mah
        weights <- wgtFUN(mah)
      }
    }
  }
  ans$mcd.wt <- rep.int(NA, length(ok))
  ans$mcd.wt[ok] <- weights
  if (names) {
    if (length(dimn[[1]])) 
      names(ans$mcd.wt) <- dimn[[1]]
    if (length(dimn[[1]])) 
      dimnames(x)[[1]] <- names(ans$mcd.wt)[ok]
    else dimnames(x) <- list(seq(along = ok)[ok], NULL)
  }
  if(fixedCenter){ ans$X <- savedx # PR
  } else { ans$X <- x }
  ans$wt <- NULL
  if (trace) 
    cat(ans$method, "\n")
  ans$raw.cnp2 <- raw.cnp2
  ans$cnp2 <- cnp2
  if (nsamp == "deterministic") 
    ans <- c(ans, mcd[c("iBest", "n.csteps", 
                        if (save.hsets) "initHsets")])
  class(ans) <- "mcd"
  if (is.list(ans$singularity)) 
    warning(paste(strwrap(.MCDsingularityMsg(ans$singularity, 
                                             ans$n.obs)), collapse = "\n"), domain = NA)
  if(fixedCenter) ans$center <- givenCenter # PR
  ans
} 




## Functions copied from robustbase

.fastmcd       <- utils::getFromNamespace(".fastmcd", "robustbase")
.MCDcnp2       <- utils::getFromNamespace(".MCDcnp2", "robustbase")
.MCDcnp2.rew   <- utils::getFromNamespace(".MCDcnp2.rew", "robustbase")
.detmcd        <- utils::getFromNamespace(".detmcd", "robustbase")
doScale        <- utils::getFromNamespace("doScale", "robustbase")
smoothWgt      <- utils::getFromNamespace("smoothWgt", "robustbase")
.MCDadaptWgt.c <- utils::getFromNamespace(".MCDadaptWgt.c", "robustbase")
.wgtFUN.covMcd <- utils::getFromNamespace(".wgtFUN.covMcd", "robustbase")
rrcov.control  <- utils::getFromNamespace("rrcov.control", "robustbase")
pasteK         <- utils::getFromNamespace("pasteK", "robustbase")
robScalefn     <- utils::getFromNamespace("robScalefn", "robustbase")
.MCDsingularityMsg <- utils::getFromNamespace(".MCDsingularityMsg", "robustbase")
getDefCtrl <- utils::getFromNamespace("getDefCtrl", "robustbase")
.MCDcons <-  utils::getFromNamespace(".MCDcons", "robustbase")
## End functions copied from robustbase

