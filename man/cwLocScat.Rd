\name{cwLocScat}
\alias{cwLocScat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate location and scatter of data with cellwise weights
}
\description{
  Computes different estimators of multivariate location
  and scatter for cellwise weighted data.
  }
\usage{
cwLocScat(X, W, methods = "all", lmin = 1e-3,
                     crit = 1e-12, maxiter= 1000, 
                     initCwCov = FALSE, initEst = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{X}{An \eqn{n} by \eqn{d} data matrix or data frame. Must be
              given. \code{X} is allowed to contain \code{NA}'s.}
  \item{W}{An \eqn{n} by \eqn{d} matrix of nonnegative cellwise weights.
              Must be given. \code{W} is not allowed to contain \code{NA}'s.}
  \item{methods}{either \code{"all"} or \code{"explicit"}. If \code{"explicit"}
              only the explicit estimates cwMean, cwCov
              and sqrtCov are computed. If \code{"all"} (the
              default) also the cellwise MLE is carried out,
              yielding cwMLEmu and cwMLEsigma.}
  \item{lmin}{if not \code{NULL}, a lower bound on the eigenvalues
              of the estimated covariance matrices on the
              standardized data, to avoid singularity.}
  \item{crit}{convergence criterion of successive mu
              and Sigma estimates in the EM algorithm.}
  \item{maxiter}{maximal number of iteration steps in EM.}
  \item{initCwCov}{if \code{TRUE}, uses the weighted mean and cwCov
              as initial estimates for the weighted EM.}
  \item{initEst}{if not \code{NULL}, a list with initial estimates \code{$mu}
              of the mean, \code{$Sigma} of the covariance matrix,
              for the weighted EM. Has no effect when
              \code{initCwCov = TRUE}.}
}

\value{

 A list with components: \cr
  \itemize{
  \item{\code{cwMean}\cr the explicit cellwise weighted mean.}
  \item{\code{cwCov}\cr explicit cellwise weighted covariance matrix.
               Is asymptotically normal but not necessarily
               PSD (unless a nonnegative \code{lmin} was specified).}
  \item{\code{sqrtCov}\cr the cellwise weighted covariance matrix of Van
               Aelst et al (2011). Also asymptotically normal
               but not necessarily PSD (unless a nonnegative
               \code{lmin} was specified).}
  \item{\code{cwMLEmu}\cr the location estimate obtained by the cwMLE.}
  \item{\code{cwMLEsigma}\cr the covariance matrix obtained by the cwMLE.
               Is PSD when the EM algorithm converges.}
  }
}


\references{
P.J. Rousseeuw (2022). Analyzing cellwise weighted data, ArXiv:2209.12697.
\href{https://arxiv.org/abs/2209.12697}{(link to open access pdf)}
}


\author{
P.J. Rousseeuw
}


\seealso{ 
\code{\link{weightedEM}},
\code{\link{unpack}}
}

\examples{
data("data_personality_traits")
X <- data_personality_traits$X
W <- data_personality_traits$W
fit <- cwLocScat(X, W)
fit$cwMLEiter # number of iteration steps taken
round(fit$cwMLEmu, 2)
round(fit$cwMean, 2)
round(fit$cwMLEsigma, 2)
round(fit$cwCov, 2)


# For more examples, we refer to the vignette:
\dontrun{
vignette("cellwise_weights_examples")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
