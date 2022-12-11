\name{data_clothes}
\alias{data_clothes}
\docType{data}
\title{
The clothes dataset
}
\description{
The clothes dataset contains a contingency table of trade flows from outside the European Union into each of its 28 member states. The columns in the contingency table in Riani et al. (2022) are five different price brackets, from lowest to highest.
}
\usage{data("data_clothes")}
\format{
  A matrix with 28 observations of 5 price brackets.
}
\source{
Riani, M., Atkinson, A. C., Torti, F., Corbellini, A. (2022). Robust correspondence analysis. \emph{Journal of the Royal Statistical Society Series C: Applied Statistics}, \bold{71}(5), 1381--1401.
}
\references{
Raymaekers and Rousseeuw (2022), Challenges of cellwise outliers.
}
\examples{
data(data_clothes)
}
