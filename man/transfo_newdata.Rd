\name{transfo_newdata}
\alias{transfo_newdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform variables based on the output of \code{\link{transfo}}. 
}
\description{
Based on the output of \code{\link{transfo}}, transform the variables using Yeo-Johnson and/or Box-Cox transformations with the previously estimated parameters and standardization.
}
\usage{
transfo_newdata(Xnew, transfo.out)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xnew}{A data matrix with d columns, which contain the variables to be transformed. The number of columns and their names must be the same as those of the original data on which \code{\link{transfo}} was run. The number of rows may be different.

}
\item{transfo.out}{
The output of a call to \code{\link{transfo}}.
}
}
\value{
Returns a matrix with transformed variables.
}
\references{
J. Raymaekers and P.J. Rousseeuw (2021). Transforming variables to central normality. \emph{Machine Learning}. \doi{10.1007/s10994-021-05960-5}(link to open access pdf)
}
\author{
J. Raymaekers and P.J. Rousseeuw
}

\seealso{
\code{\link{transfo}}
}

\examples{
set.seed(123); tempraw <- matrix(rnorm(2000), ncol = 2)
tempx <- cbind(tempraw[, 1],exp(tempraw[, 2]))
tempy <- 0.5 * tempraw[, 1] + 0.5 * tempraw[, 2] + 1
x <- tempx[1:900, ]
y <- tempy[1:900]
tx.out <- transfo(x, type = "bestObj")
tx.out$ttypes
tx.out$lambdahats
tx <- tx.out$Y
lm.out <- lm(y ~ tx)
summary(lm.out)
xnew <- tempx[901:1000, ]
xtnew <- transfo_newdata(xnew, tx.out)
yhatnew <- tcrossprod(lm.out$coefficients, cbind(1, xtnew)) 
plot(tempy[901:1000], yhatnew); abline(0, 1)
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
