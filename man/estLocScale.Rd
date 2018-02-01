\name{estLocScale}
\alias{estLocScale}

\title{
Estimate robust location and scale
}
\description{
Estimate a robust location estimate and scale estimate of every column in \code{X}.
}
\usage{
estLocScale(X, type = "wrap", precScale = 1e-12)
}

\arguments{
  \item{X}{
  The input data. It must be an \eqn{n} by \eqn{d} matrix or a data frame.
}
  \item{type}{
 The type of estimators used. One of:
 \itemize{
  \item{"mcd": \cr the location is the weighted univariate MCD estimator with cutoff sqrt(qchisq(0.975,1)). The scale is the corresponding weighted univariate MCD estimator, with a correction factor to make it approximately unbiased at gaussian data.}
 \item{"wrap": \cr Starting from the initial estimates corresponding to option \code{"mcd"}, the location is the 1-step M-estimator with the wrapping psi function with \eqn{b = 1.5} and \eqn{c = 4}. The scale estimator is the same as in option \code{"mcd"}.
 }

 }
 Defaults to "wrap".
}

  \item{precScale}{
The precision scale used throughout the algorithm. Defaults to \eqn{1e-12}.
}
}


\value{
A list with components: \cr
  \itemize{
    \item{\code{loc} \cr
   A vector with the estimated locations.
    }
    \item{\code{scale} \cr
   A vector with the estimated scales.
    }
}
}

\references{
 Raymaekers, J., Rousseeuw P.J. (2018). Fast robust correlation for high dimensional data.
 \emph{arXiv:1712.05151}
}
\author{
 Raymaekers, J. and Rousseeuw P.J.
}

\seealso{
\code{\link{wrap}}
}

\examples{
library(MASS) 
set.seed(12345) 
n = 100; d = 10
X = mvrnorm(n, rep(0, 10), diag(10))
locScale = estLocScale(X)
}