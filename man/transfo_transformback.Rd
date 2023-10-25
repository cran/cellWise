\name{transfo_transformback}
\alias{transfo_transformback}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Backtransform variables based on the output of \code{\link{transfo}}. 
}
\description{
Based on the output of \code{\link{transfo}}, backtransform the variables to their original shape through the inverse Yeo-Johnson and/or Box-Cox transformations with the previusly estimated parameters and standardization.
}
\usage{
transfo_transformback(Ynew, transfo.out)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ynew}{A data matrix with d columns, which contain the variables to be backtransformed. The number of columns must be the same as the output \code{Y} of the run of \code{\link{transfo}} on the  original data. The number of rows may be different.

}
\item{transfo.out}{
The output of a call to \code{\link{transfo}}.
}
}

\value{
Returns a matrix with backtransformed variables.
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
set.seed(123); x <- matrix(rnorm(2000), ncol = 2)
y <- sqrt(abs(0.3 * x[, 1] + 0.5 * x[, 2] + 4))
ty.out <- transfo(y, type = "BC")
ty.out$lambdahats
ty <- ty.out$Y
lm.out <- lm(ty ~ x)
yhat <- transfo_transformback(lm.out$fitted.values, ty.out)
plot(y, yhat); abline(0, 1)
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
