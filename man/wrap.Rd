\name{wrap}
\alias{wrap}
\title{
Wrap the data.
}
\description{
Transforms multivariate data X using the wrapping function with b = 1.5 and c = 4 and the location and scale given in locX and scaleX. 
}
\usage{
wrap(X, locX, scaleX, precScale = 1e-12)
}

\arguments{
  \item{X}{the input data. It must be an \eqn{n} by \eqn{d} matrix or a data frame.
}
  \item{locX}{The location estimates of the columns of the input data \code{X}. Must be a vector of length \eqn{d}.
}
  \item{scaleX}{The scale estimates of the columns of the input data \code{X}. Must be a vector of length \eqn{d}.
}
  \item{precScale}{The precision scale used throughout the algorithm. Defaults to \eqn{1e-12}
}
}

\value{A list with components: \cr
  \itemize{
    \item{\code{Xw} \cr
   The wrapped data.
    }
    \item{\code{colInWrap} \cr 
    The column numbers for which the scale estimate was larger than precScale. Those with scale estimate <= precScale do not occur in Xw to avoid division by (near) zero.
    }
}
}
\references{
  Raymaekers, J., Rousseeuw P.J. (2019). Fast robust correlation for high dimensional data. \emph{Technometrics}, published online.
}
\author{
 Raymaekers, J. and Rousseeuw P.J.
}

\seealso{
\code{\link{estLocScale}}
}
\examples{
library(MASS) 
set.seed(12345) 
n <- 100; d <- 10
X <- mvrnorm(n, rep(0, 10), diag(10))
locScale <- estLocScale(X)
Xw <- wrap(X, locScale$loc, locScale$scale)$Xw
}
