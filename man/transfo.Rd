\name{transfo}
\alias{transfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robustly fit the Box-Cox or Yeo-Johnson transformation 
}
\description{
This function uses reweighted maximum likelihood to robustly fit the 
Box-Cox or Yeo-Johnson transformation to each variable in a dataset. 
Note that this function first calls \code{\link{checkDataSet}} to ensure that the variables to be transformed are not too discrete.
}
\usage{
transfo(X, type = "YJ", robust = TRUE,
        standardize = TRUE,
        quant = 0.99, nbsteps = 2, checkPars = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ A data matrix of dimensions n x d.
  Its columns are the variables to be transformed.
}
\item{type}{ The type of transformation to be fit. Should be one of: \itemize{
\item{\code{"BC"}}: Box-Cox power transformation. Only works
 for strictly positive variables. If this type
 is given but a variable is not strictly
  positive, the function stops with a
   message about that variable.
\item{\code{"YJ"}} Yeo-Johnson power transformation. The data
 may have positive as well as negative values.
\item{\code{"bestObj"}} for strictly positive variables both BC
                         and YJ are run, and the solution with
   lowest objective is kept. On the other
  variables YJ is run.
}

}
\item{robust}{ if \code{TRUE} the Reweighted Maximum 
                Likelihood method is used, which first
   computes a robust initial estimate of the
 transformation parameter lambda. If \code{FALSE} the classical ML method is used.
 }
 \item{standardize}{ whether to standardize the variables \strong{before and after} the
 power transformation. See Details below.
 }
 \item{quant}{quantile for determining the weights in the 
   reweighting step (ignored when \code{robust=FALSE}).
   }
 \item{nbsteps}{number of reweighting steps (ignored when 
   \code{robust=FALSE}).
   }
 \item{checkPars}{Optional list of parameters used in the call to
 \code{\link{checkDataSet}}. The options are:
  \itemize{
  \item{ \code{coreOnly} \cr
      If \code{TRUE}, skip the execution of checkDataset. Defaults to \code{FALSE}}
       \item{\code{numDiscrete}\cr
  A column that takes on numDiscrete or fewer values
  will be considered discrete and not retained in the cleaned data.
  Defaults to \eqn{5}.
}
  \item{\code{precScale} \cr
  Only consider columns whose scale is larger than precScale.
  Here scale is measured by the median absolute deviation.
  Defaults to \eqn{1e-12}.
}
 \item{\code{silent}\cr
  Whether or not the function progress messages should be printed.
  Defaults to \code{FALSE}.
}
      }
   }
}
\value{
A list with components: \cr
  \itemize{
    \item{\code{lambdahats} \cr
    the estimated transformation parameter for each column of \code{X}.
    }
    \item{\code{Y} \cr
    A matrix in which each column is the transformed version of the
   corresponding column of \code{X}.
   The transformed version includes pre- and post-standardization if \code{standardize=TRUE}.
    }
      \item{\code{muhat} \cr
    The estimated location of each column of \code{Y}.
    }
      \item{\code{sigmahat} \cr
    The estimated scale of each column of \code{Y}.
    }
     
     \item{\code{weights} \cr
    The final weights from the reweighting.
    }
      \item{\code{ttypes} \cr
    The type of transform used in each column.
    }
     \item{\code{objective} \cr
     Value of the (reweighted) maximum likelihood objective function.
    }
    \item{values of \code{\link{checkDataSet}}, unless \code{coreOnly} is \code{TRUE}.}
    }
}

\details{
In case \code{standardize = TRUE}, the variables is standardized before and after transformation. 
For BC the variable is divided by its median before transformation.
For YJ and \code{robust = TRUE} this subtracts its median and divides by its mad (median absolute deviation) before transformation. For YJ and \code{robust = FALSE} this subtracts the mean and divides by the standard deviation before transformation. For the standardization after the transformation, the classical mean and standard deviation are used in case  \code{robust = FALSE}. If \code{robust = TRUE}, the mean and standard deviation are calculated robustly on a subset of inliers.  
}
\references{
J. Raymaekers and P.J. Rousseeuw (2021). Transforming variables to central normality. \emph{Machine Learning}. \doi{10.1007/s10994-021-05960-5}(link to open access pdf)
}
\author{
J. Raymaekers and P.J. Rousseeuw
}
\seealso{
\code{\link{transfo_newdata}}, \code{\link{transfo_transformback}}
}

\examples{

# find Box-Cox transformation parameter for lognormal data:
set.seed(123)
x <- exp(rnorm(1000))
transfo.out <- transfo(x, type = "BC")
# estimated parameter:
transfo.out$lambdahat
# value of the objective function:
transfo.out$objective
# the transformed variable:
transfo.out$Y
# the type of transformation used:
transfo.out$ttypes
# qqplot of the transformed variable:
qqnorm(transfo.out$Y); abline(0,1)

# For more examples, we refer to the vignette:
\dontrun{
vignette("transfo_examples")
}
}