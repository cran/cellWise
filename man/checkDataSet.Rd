\name{checkDataSet}
\alias{checkDataSet}

\title{
Purify the dataset
}

\description{
This function checks the dataset X, and sets aside certain
columns and rows that do not satisfy the conditions.
It is used by the \code{\link{DetectDeviatingCells}} function but can be used by itself, to clean a dataset for a different type of analysis.
}
\usage{
checkDataSet(X, fracNA = 0.5, numDiscrete = 3, precScale = 1e-12)
}

\arguments{
  \item{X}{\code{X} is the input data, and must be an \eqn{n} by \eqn{d} matrix or data frame.
}
  \item{fracNA}{
   Only retain columns and rows with fewer NAs than this fraction.
   Defaults to \eqn{0.5}.
}
  \item{numDiscrete}{
  A column that takes on numDiscrete or fewer values
  will be considered discrete and not retained in the cleaned data.
  Defaults to \eqn{3}.
}
  \item{precScale}{
  Only consider columns whose scale is larger than precScale.
  Here scale is measured by the median absolute deviation.
  Defaults to \eqn{1e-12}.
}
}

\value{
A list with components: \cr
  \itemize{
    \item{\code{colInAnalysis} \cr
    Column indices of the columns used in the analysis.
    }
    \item{\code{rowInAnalysis}} \cr
    Row indices of the rows used in the analysis.
    \item{\code{namesNotNumeric} \cr
    Names of the variables which are not numeric.       
    }
    \item{\code{namesCaseNumber} \cr
    Names of the cases.
    }
    \item{\code{namesNAcol} \cr
    Names of the NA columns.          
    }
    \item{\code{namesNArow} \cr
    Names of the NA rows.           
    }
    \item{\code{namesDiscrete} \cr
    Names of the discrete variables.           
    }
    \item{\code{namesZeroScale} \cr
    Names of the variables with zero scale.           
    }
    \item{\code{remX} \cr
    Cleaned data after checkDataSet.           
    }
  }
}
\references{
 Rousseeuw P.J., Van den Bossche W. (2016). Detecting deviating data cells.
 \emph{arXiv:1601.07251}
}

\author{
 Rousseeuw P.J., Van den Bossche W. 
}

\seealso{
\code{\link{DetectDeviatingCells}}
}

\examples{
library(MASS) 
set.seed(12345) 
n = 100; d = 10
A = matrix(0.9, d, d); diag(A) = 1
xclean = mvrnorm(n, rep(0,d), A)
checkedx = checkDataSet(xclean)

# For more examples, we refer to the vignette:
vignette("DDC_examples")
}
