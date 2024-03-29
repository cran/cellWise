\name{MacroPCApredict}

\alias{MacroPCApredict}

\title{
MacroPCApredict
}

\description{Based on a \code{\link{MacroPCA}} fit of an initial (training) data set \code{X}, this function analyzes a
new (test) data set \code{Xnew}.}

\usage{
MacroPCApredict(Xnew, InitialMacroPCA, MacroPCApars = NULL)
}

\arguments{

\item{Xnew}{The new data (test data), which must be a matrix or a data frame.
It must always be provided. Its columns (variables) should correspond to those of \code{InitialMacroPCA$remX}.}

\item{InitialMacroPCA}{The output of the MacroPCA function on the initial
(training) dataset. Must be provided.}

\item{MacroPCApars}{The input options to be used for the prediction.
By default the options of InitialMacroPCA are used. For the complete list of
options see the function \code{\link{MacroPCA}}.}
  
}

\value{
  A list with components: \cr
  \item{MacroPCApars}{the options used in the call.}
  \item{DDC}{result of DDCpredict which is the first step of MacroPCApredict.
  See the function \code{\link{DDCpredict}}.}
  \item{scaleX}{the scales of the columns of \code{X}.}
  \item{k}{the number of principal components.}
  \item{loadings}{the columns are the \code{k} loading vectors.}
  \item{eigenvalues}{the \code{k} eigenvalues.}
  \item{center}{vector with the fitted center.}
  \item{It}{number of iteration steps.}
  \item{diff}{convergence criterion.}
  \item{Xnew.NAimp}{\code{Xnew} with all \code{NA}'s imputed by \code{MacroPCA}.}
  \item{scores}{scores of \code{Xnew.NAimp}.}
  \item{OD}{orthogonal distances of the rows of \code{Xnew.NAimp}.}
  \item{cutoffOD}{cutoff value for the OD.}
  \item{SD}{score distances of the rows of \code{Xnew.NAimp}.}
  \item{cutoffSD}{cutoff value for the SD.}
  \item{highOD}{row numbers of cases in \code{Xnew.NAimp} whose \code{OD} is above \code{cutoffOD}.}
  \item{highSD}{row numbers of cases in \code{Xnew.NAimp} whose \code{SD} is above \code{cutoffSD}.}
  \item{residScale}{scale of the residuals.}
  \item{stdResid}{standardized residuals. Note that these are \code{NA} for all missing values of \code{Xnew}.}
  \item{indcells}{indices of cellwise outliers.}
  \item{NAimp}{various results for the NA-imputed Xnew.}
  \item{Cellimp}{various results for the cell-imputed Xnew.}
  \item{Fullimp}{various result for the fully imputed Xnew.}
}

\references{
Hubert, M., Rousseeuw, P.J., Van den Bossche W. (2019). MacroPCA: An all-in-one PCA method allowing for missing values as well as cellwise and rowwise outliers. \emph{Technometrics}, \bold{61}(4), 459-473. \href{https://www.tandfonline.com/doi/full/10.1080/00401706.2018.1562989}{(link to open access pdf)}
}

\author{
 Rousseeuw P.J., Van den Bossche W. 
}


\seealso{
\code{\link{checkDataSet}}, \code{\link{cellMap}},
\code{\link{DDC}}, \code{\link{DDCpredict}},
\code{\link{MacroPCA}}
}

\examples{
library(MASS) 
set.seed(12345) 
n <- 50; d <- 10
A <- matrix(0.9, d, d); diag(A) = 1
x <- mvrnorm(n, rep(0,d), A)
x[sample(1:(n * d), 50, FALSE)] <- NA
x[sample(1:(n * d), 50, FALSE)] <- 10
MacroPCA.out <- MacroPCA(x, 2)
xnew <- mvrnorm(25, rep(0,d), A)
xnew[sample(1:(25 * d), 12, FALSE)] <- 10
predict.out <- MacroPCApredict(xnew, MacroPCA.out)
cellMap(predict.out$stdResid)

# For more examples, we refer to the vignette:
\dontrun{
vignette("MacroPCA_examples")
}
}
