\name{plot_cellMCD}
\alias{plot_cellMCD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw plots based on the cellwise minimum covariance determinant estimator cellMCD
}
\description{
Function for making plots based on the output of \code{\link{cellMCD}}.
}
\usage{
plot_cellMCD(cellout, type = "Zres/X", whichvar = NULL,
             horizvar = NULL, vertivar = NULL,  
             hband = NULL, vband = NULL, drawellipse = T,
             opacity = 0.5, identify = FALSE, 
             ids = NULL, labelpoints = T, vlines = FALSE,
             clines = TRUE, main = NULL,
             xlab = NULL, ylab = NULL, xlim = NULL,
             ylim = NULL, cex = 1, cex.main = 1.2, 
             cex.txt = 0.8, cex.lab = 1, line = 2.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cellout}{output of function \code{\link{cellMCD}}}
  \item{type}{type of diagnostic plot. Should be one of \code{"index"}, \code{"Zres/X"}, \code{"Zres/pred"}, \code{"X/pred"}, or \code{"bivariate"}.}
  \item{whichvar}{ number or name of the variable to be plotted. Not applicable when \code{type = "bivariate"}.}
  \item{horizvar}{ number or name of the variable to be plotted on the horizontal axis. Only when \code{type = "bivariate"}.}
  \item{vertivar}{ number or name of the variable to be plotted on the  vertical axis. Only when \code{type = "bivariate"}.}
  \item{hband}{ whether to draw a horizontal tolerance band. \code{TRUE} or \code{FALSE}. \code{NULL} yields \code{TRUE} when \code{type} is \code{"index"}, \code{"Zres/X"}, or \code{"Zres/pred"}.}
  
  \item{vband}{ whether to draw a vertical tolerance band. \code{TRUE} or \code{FALSE}.  \code{NULL} yields \code{TRUE} when \code{type} is \code{"Zres/X"} or \code{"Zres/pred"}.}
  \item{drawellipse}{ whether to draw a \eqn{99}\% tolerance ellipse. Only for \code{type = "bivariate"}.}
  \item{opacity}{opacity of the plotted points: 1 is fully opaque, less is more transparent.}
  \item{identify}{ if \code{TRUE}, identify cases by mouseclick, then \code{Esc}.}
  \item{ids}{vector of case numbers to be emphasized (colored red) in the plot. If \code{NULL} or of length zero, none are emphasized.}
  \item{labelpoints}{ if \code{TRUE}, labels the points in ids by their row name in \code{X}.}
  \item{vlines}{ for the points in \code{ids}, draw dashed vertical lines from their standardized residual to 0 when \code{type} is \code{"index"}, \code{"Zres/X"}, or \code{"Zres/pred"}. Draws dashed vertical lines to the diagonal when \code{type = "X/pred"}. Can be \code{TRUE} or \code{FALSE}, default is \code{FALSE}.}
  \item{clines}{ only for type == "bivariate". If TRUE, draws
  a red connecting line from each point in ids to 
  its imputed point, shown in blue.}
  
  
  \item{main}{     main title of the plot. If \code{NULL}, it is constructed automatically from the arguments.}  
  \item{xlab}{ overriding label for x-axis, unless \code{NULL}.}
  \item{ylab}{ overriding label for y-axis, unless \code{NULL}.}
  \item{xlim}{ overriding limits of horizontal axis.}
  \item{ylim}{ overriding limits of vertical axis.}
  \item{cex}{  size of plotted points.}
  \item{cex.main}{size of the main title.}
  \item{cex.lab}{ size of the axis labels.}
  \item{cex.txt}{ size of the point labels.}
  \item{line}{ distance of axis labels to their axis.}
}
   
  
  
              
\value{
 \code{NULL}, unless \code{identify = TRUE}. Then a list with components: \cr
  \itemize{
   \item{\code{ids} \cr
    the case number(s) that were identified
    }
    \item{\code{coords} \cr
    coordinates of all points in the plot.
    }
    }
}
\references{
J. Raymaekers and P.J. Rousseeuw (2022). The cellwise MCD estimator, Journal of the American Statistical Association, to appear.
\doi{10.1080/01621459.2023.2267777}(link to open access pdf)
}
\author{
J. Raymaekers and P.J. Rousseeuw
}


\seealso{
\code{\link{cellMCD}}
}
\examples{
mu <- rep(0, 3)
Sigma <- diag(3) * 0.5 + 0.5
set.seed(123)
X <- MASS::mvrnorm(1000, mu, Sigma)
X[1:5, 1] <- X[1:5, 1] + 5
X[6:10, 2] <- X[6:10, 2] - 10
X[12, 1:2] <- c(-4,8)
cellMCD.out <- cellMCD(X)
plot_cellMCD(cellMCD.out, type="bivariate", 
             horizvar=1, vertivar=2, ids=c(1:10,12))

# For more examples, we refer to the vignette:
\dontrun{
vignette("cellMCD_examples")
}
}