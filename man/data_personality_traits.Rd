\name{data_personality_traits}
\alias{data_personality_traits}
\docType{data}
\title{
The personality traits data
}
\description{
This dataset describes personality traits of 10 persons. The variables are the 6 traits Anxiety, Agoraphobia, Arachnophobia, Adventurous, Extraversion, and Sociability. 
}
\usage{data("data_personality_traits")}
\format{
The data contains a list with two elements:
\itemize{
\item{\code{X}\cr
a \eqn{10} by \eqn{6} matrix of values describing \eqn{6} personality traits for each of the \eqn{10} participants. 
}
\item{\code{W} \cr
a \eqn{10} by \eqn{6} matrix of cellwise weights. Each weight is the inverse of the length of the support of the membership function of the fuzzy number in the original data set.
}
}
}
\source{
G. Hesamian, and Akbari, M. G. (2019), 
Principal component analysis based on intuitionistic fuzzy random variables, \emph{Computational and Applied Mathematics}, \bold{38}(158), 1--14.
}
\references{
P.J. Rousseeuw (2022). Analyzing cellwise weighted data, ArXiv:2209.12697.
\href{https://arxiv.org/abs/2209.12697}{(link to open access pdf)}
}



\examples{
data(data_personality_traits)

# For the examples in Rousseeuw (2022), see:
\dontrun{
vignette("cellwise_weights_examples")
}
}
