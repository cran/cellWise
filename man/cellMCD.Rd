\name{cellMCD}
\alias{cellMCD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cellWise minimum covariance determinant estimator
}
\description{
The cellwise minimum covariance determinant estimator
computes cellwise robust estimates of the center and covariance matrix of a data set \code{X}. The algorithm guarantees a monotone decrease of an objective function, 
which is based on observed Gaussian log-likelihood. By default, it starts by calling \code{\link{checkDataSet}} to clean the data.
}
\usage{
cellMCD(X, alpha = 0.75, quant = 0.99,
        crit = 1e-4, noCits = 100, lmin = 1e-4,
        checkPars = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{\code{X} is the input data, and must be an \eqn{n} by \eqn{d} matrix or a data frame.
  }
  \item{alpha}{In each column, at least \eqn{n*}\code{alpha} cells must remain unflagged. Defaults to \eqn{75}\%, should not be set (much) lower.
  }
  \item{quant}{Determines the cutoff value to flag cells. Defaults to \eqn{0.99}.
  }
  
  \item{crit}{The iteration stops when successive covariance matrices (of the standardized data) differ by less than \code{crit}. Defaults to \eqn{1e-4}.
}

\item{noCits}{The maximal number of C-steps used.}

 \item{lmin}{ a lower bound on the eigenvalues of the estimated covariance matrix on the standardized data. Defaults to \eqn{1e-4}. Should not be smaller than \eqn{1e-6}.}

  \item{checkPars}{Optional list of parameters used in the call to
 \code{\link{checkDataSet}}. The options are:
  \itemize{
  \item{ \code{coreOnly} \cr
      If \code{TRUE}, skip the execution of checkDataset. Defaults to \code{FALSE}.}
       \item{\code{numDiscrete}\cr
  A column that takes on numDiscrete or fewer values
  will be considered discrete and not retained in the cleaned data.
  Defaults to \eqn{5}.
}
\item{fracNA}{\cr
   Only retain columns and rows with fewer NAs than this fraction.
   Defaults to \eqn{0.5}.
}
  \item{\code{precScale} \cr
  Only consider columns whose scale is larger than precScale.
  Here scale is measured by the median absolute deviation.
  Defaults to \eqn{1e-12}.
}
 \item{\code{silent}\cr
  Whether or not the function progress messages should be suppressed.
  Defaults to \code{FALSE}.
}
      }
   }
  }
  
              
\value{
 A list with components: \cr
  \itemize{
    \item{\code{mu} \cr
    the cellMCD estimate of location.
    }
    \item{\code{S} \cr
    the cellMCD estimate of scatter.
    }
    \item{\code{W} \cr
    the cellMCD estimate of \code{W}, a binary matrix indicating all outlying cells as zero.
    }
    \item{\code{preds} \cr
    predictions (=conditional expectations) of the flagged cells, given the clean cells in the same row.
    }
    \item{\code{csds} \cr
    conditional standard deviations of the flagged cells, given the clean cells in the same row.
    }
    \item{\code{Ximp} \cr
    imputed data matrix.
    }
    \item{\code{Zres} \cr
    matrix of cellwise standardized residuals.
    }
    \item{\code{rscales} \cr
    robust scales used to standardize the data before running the algorithm. The results \code{m}, \code{S}, \code{preds}, \code{Ximp} are returned in their original scale. 
    }
    \item{\code{nosteps} \cr
    number of steps the algorithm took to converge.
    }
     \item{\code{X} \cr
    the data on which the algorithm was executed.
    }
    \item{\code{quant} \cr
    the cutoff used to flag the cells.
    }
    }
}
\references{
J. Raymaekers and P.J. Rousseeuw (2022). The cellwise MCD estimator, ArXiv:2207.13493.
\href{https://arxiv.org/abs/2207.13493}{(link to open access pdf)}
}
\author{
J. Raymaekers and P.J. Rousseeuw
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot_cellMCD}}
}
\examples{
mu    <- rep(0, 3)
Sigma <- diag(3) * 0.5 + 0.5
set.seed(123)
X <- MASS::mvrnorm(1000, mu, Sigma)
X[1:5, 1]  <- X[1:5, 1] + 5
X[6:10, 2] <- X[6:10, 2] - 10
X[12, 1:2] <- c(-4,8)
cellMCD.out <- cellMCD(X)
cellMCD.out$mu
cov2cor(cellMCD.out$S)
cellMCD.out$W[1:15,]
cellMCD.out$Ximp[1:15,]
cellMap(D=X[1:15,], R=cellMCD.out$Zres[1:15,],
        columnlabels=c("X1","X2","X3"),
        rowlabels=1:15, drawCircles=FALSE)
# For more examples, we refer to the vignette:
\dontrun{
vignette("cellMCD_examples")
}
}