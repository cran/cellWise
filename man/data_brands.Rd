\name{data_brands}
\alias{data_brands}
\docType{data}
\title{
The brands dataset
}
\description{
The brands data is a contingency table summarizing the 2014 Auto Brand Perception survey by Consumer Reports (USA), which is publicly available on https://boraberan.wordpress.com/2016/09/22/. The survey questioned 1578 participants on what they considered attributes of 39 different car brands.
}
\usage{data("data_brands")}
\format{
  A matrix with 39 observations of 7 attributes. The attributes (columns) are Fuel Economy, Innovation, Performance, Quality, Safety, Style and Value.
}
\source{
https://boraberan.wordpress.com/2016/09/22/.
}
\references{
Riani, M., Atkinson, A. C., Torti, F., Corbellini, A. (2022). Robust correspondence analysis. \emph{Journal of the Royal Statistical Society Series C: Applied Statistics}, \bold{71}(5), 1381--1401.

Raymaekers and Rousseeuw (2022), Challenges of cellwise outliers.
}
\examples{
data(data_brands)
}
