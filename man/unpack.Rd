\name{unpack}
\alias{unpack}

\title{
Unpacks cellwise weighted data
}
\description{
This function transforms a dataset X with cellwise weights W
to an extended data matrix U with the same number of columns but more rows, and containing more NA's. Its rows have the case weights v.
}
\usage{
unpack(X,W)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{An \eqn{n} by \eqn{d} data matrix or data frame. Must be given.
      \code{X} is allowed to contain \code{NA}'s.}
  \item{W}{An \eqn{n} by \eqn{d} matrix of nonnegative cellwise weights.
      Must be given. \code{W} is not allowed to contain \code{NA}'s.}
}

\value{
 A list with components: \cr
  \itemize{
  \item{\code{U}\cr
  unpacked data matrix, with the same columns as \code{X} but typically more rows.}
  \item{\code{V}\cr vector with the rowwise (=casewise) weights of \code{U}.}
  }
}

\references{
P.J. Rousseeuw (2022). Analyzing cellwise weighted data, ArXiv:2209.12697.
\href{https://arxiv.org/abs/2209.12697}{(link to open access pdf)}
}


\author{
P.J. Rousseeuw
}


\seealso{ 
\code{\link{weightedEM}},
\code{\link{cwLocScat}}
}

\examples{


X <- matrix(c(2.8, 5.3, 4.9, 7.4,
             2.3, 5.7, 4.3, 7.2,
             2.5, 5.1, 4.4, 7.6), nrow = 3, byrow = TRUE)
W <- matrix(c(0.8, 1.0, 0.3, 0.4, 
             0.3, 0.5, 0.9, 0.5, 
             1.0, 0.6, 0, 0.7), nrow = 3, byrow = TRUE)
rownames(X) <- rownames(W) <- c("A", "B", "C")
colnames(X) <- colnames(W) <- c("V1", "V2", "V3", "V4")
X
W
out <- unpack(X, W)
cbind(out$U, out$v)


# For more examples, we refer to the vignette:
\dontrun{
vignette("cellwise_weights_examples")
}
}