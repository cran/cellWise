\name{weightedEM}
\alias{weightedEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates location and scatter on incomplete data with case weights
}
\description{
Carries out a rowwise weighted EM algorithm to estimate mu and Sigma of incomplete Gaussian data.
}
\usage{
weightedEM(X, w=NULL, lmin=NULL, crit=1e-4, 
                      maxiter=1000, initEst=NULL, computeloglik=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{n by d data matrix or data frame.}
  \item{w}{vector with n nonnegative rowwise (casewise)
weights. If \code{NULL}, all weights are set to 1 so
an unweighted EM is carried out.}
  \item{lmin}{if not \code{NULL}, a lower bound on the eigenvalues
of the estimated EM covariance matrix on the
standardized data, to avoid singularity.}
  \item{crit}{convergence criterion of successive mu
and Sigma estimates.}
  \item{maxiter}{maximal number of iteration steps.}
  \item{initEst}{if not \code{NULL}, a list with initial estimates \code{$mu} of the mean, \code{$Sigma} of the covariance matrix.}
  \item{computeloglik}{if \code{TRUE}, the log(likelihood) is computed in
every step and reported. Default is \code{FALSE}
to save computation time.}
  
}
\value{

 A list with components: \cr
  \itemize{
  \item{\code{mu}\cr the estimated location vector.}
  \item{\code{Sigma}\cr the estimated covariance matrix.}
  \item{\code{impX}\cr the imputed data matrix.}
  \item{\code{niter}\cr the number of iteration steps taken.}
  \item{\code{loglikhd}\cr vector with the total log(likelihood) at every
             iteration step. When \code{computeloglik = FALSE} this
             array contains NA's.}
  }

}
\references{
P.J. Rousseeuw (2022). Analyzing cellwise weighted data, ArXiv:2209.12697.
\href{https://arxiv.org/abs/2209.12697}{(link to open access pdf)}
}


\author{
P.J. Rousseeuw
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ 
\code{\link{unpack}},
\code{\link{cwLocScat}}
}

\examples{

Sigma <- matrix(0.7, 3, 3); diag(Sigma) <- 1
set.seed(12345); X <- MASS::mvrnorm(1000, rep(0, 3), Sigma)
X[1, 3] <- X[2, 2] <- X[3, 1] <- X[4, 1] <- X[5, 2] <- NA
w <- runif(1000, 0, 1) # rowwise weights
out <- weightedEM(X, w, crit = 1e-12, computeloglik = TRUE)
out$niter # number of iteration steps taken
plot(1:out$niter, out$loglikhd[1:out$niter], type = 'l',
     lty = 1, col = 4, xlab = 'step', ylab = 'log(likelihood)',
     main = 'log(likelihood) of weighted EM iterations')
out$mu # estimated center
round(out$Sigma, 6) # estimated covariance matrix
head(X) # the data has NA's
head(out$impX) # imputed data, has no NA's

# For more examples, we refer to the vignette:
\dontrun{
vignette("cellwise_weights_examples")
}
}
