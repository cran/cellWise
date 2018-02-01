\name{DetectDeviatingCells}

\alias{DetectDeviatingCells}

\title{
Detect Deviating Cells
}

\description{
This function aims to detect cellwise outliers in the data. These are entries in the data matrix which are substantially higher or lower than what could be expected based on the other cells in its column as well as the other cells in its row, taking the relations between the columns into account. (Note that this function first calls \code{\link{checkDataSet}} and analyzes the remaining cleaned data.)
}

\usage{
DetectDeviatingCells(X, DDCpars = list())
}

\arguments{
  \item{X}{\code{X} is the input data, and must be an \eqn{n} by \eqn{d} matrix or a data frame.
  }
  \item{DDCpars}{ A list of available options:
    \itemize{
     \item{ \code{fastDDC} \cr
      Whether to use the fastDDC option or not. Defaults to \eqn{TRUE}.}
      \item{ \code{transFun} \cr
      The transformation function used to compute the robust correlations when \code{fastDDC} \eqn{ = TRUE}. Can be \code{"wrap"} or \code{"rank"}. Defaults to \code{"wrap"}.}
      \item{ \code{fracNA} \cr
      Only consider columns and rows with fewer NAs (missing
      values) than this fraction (percentage). Defaults to \eqn{0.5}.}
      \item{ \code{numDiscrete} \cr
      A column that takes on \code{numDiscrete} or fewer values will
      be considered discrete and not used in the analysis. Defaults to \eqn{3}.}
      \item{ \code{precScale} \cr
      Only consider columns whose scale is larger than \code{precScale}.
      Here scale is measured by the median absolute deviation. Defaults to \eqn{1e-12}.}  
      \item{ \code{tolProb} \cr
      Tolerance probability, with default \eqn{0.99}, which
      determines the cutoff values for flagging outliers.
      Used in several the steps of the algorithm.}
      \item{ \code{corrlim} \cr
      When trying to estimate \eqn{z_{ij}} from other variables \eqn{h}, we 
      will only use variables \eqn{h} with \eqn{|\rho_{j,h}| \ge corrlim}.
      Variables \eqn{j} without any correlated variables \eqn{h} satisfying 
      this are considered standalone, and treated on their own. Defaults to \eqn{0.5}.}
      \item{ \code{combinRule} \cr
      The operation to combine estimates of \eqn{z_{ij}} coming from
      other variables \eqn{h}: can be \code{wmean}, \code{wmedian}, \code{mean}, \code{median}.       Defaults to \code{wmean}.}
      \item{ \code{includeSelf} \cr
      Whether or not the combination rule will include the
      variable \eqn{j} itself. Defaults to \code{TRUE}.}
      \item{ \code{rowdetect} \cr
      Whether the rule for flagging rows is to be applied. Defaults to \code{TRUE}.}
      \item{ \code{returnBigXimp} \cr
      If TRUE, the imputed data matrix \code{Ximp} in the output
      will include the rows and columns that were not
      part of the analysis (and can still contain NAs). Defaults to \code{FALSE}.}
      \item{ \code{silent} \cr
      If TRUE, statements tracking the algorithm's progress will not be printed. Defaults to \code{FALSE}.}
    }
  }
}

\value{ A list with components: \cr
  \itemize{
    \item{\code{colInAnalysis} \cr
    The column indices of the columns used in the analysis.
    }
    \item{\code{rowInAnalysis} \cr
    The row indices of the rows used in the analysis.
    }
    \item{\code{namesNotNumeric} \cr
    The names of the variables which are not numeric.      
    }
    \item{\code{namesCaseNumber} \cr
    The names of the cases which are used in the analysis.           
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
    Cleaned data after \code{checkDataSet}.           
    }
    \item{\code{Z} \cr
    Standardized \code{X}.           
    }
    \item{\code{k} \cr
    Number of neighbors used in estimation.          
    }
    \item{\code{ngbrs} \cr
    Indicates neighbors of each column, i.e. the columns most correlated with it.   
    }
    \item{\code{robcors} \cr
    Robust correlations.           
    }
    \item{\code{robslopes} \cr
    Robust slopes.          
    }              
    \item{\code{Xest} \cr
    Predicted \code{X}.           
    }
    \item{\code{stdResid} \cr
    Residuals of orginal \code{X} minus the estimated \code{Xest},standardized by       column..           
    }
    \item{\code{indcells} \cr
    Indices of the cells which were flagged in the analysis.
    }
    \item{\code{Ti} \cr
    Outlyingness (test) value of each row.         
    }
    \item{\code{indrows} \cr
    Indices of the rows which were flagged in the analysis.
    }
    \item{\code{indall} \cr
    Indices of all cells which were flagged in the analysis plus all cells in flagged rows.
    }            
    \item{\code{Ximp} \cr
    Imputed \code{X}.           
    }
  }
}

\references{
 Rousseeuw P.J., Van den Bossche W. (2016). Detecting deviating data cells.
 \emph{arXiv:1601.07251}
  Raymaekers, J., Rousseeuw P.J. (2018). Fast robust correlation for high dimensional data.
 \emph{arXiv:1712.05151}
}

\author{
 Rousseeuw P.J., Van den Bossche W. 
}


\seealso{
\code{\link{checkDataSet}},\code{\link{cellMap}}
}

\examples{
library(MASS) 
set.seed(12345) 
n = 100; d = 10
A = matrix(0.9, d, d); diag(A) = 1
xclean = mvrnorm(n, rep(0,d), A)
DDCxclean = DetectDeviatingCells(xclean)

# For more examples, we refer to the vignette:
vignette("DDC_examples")
}
