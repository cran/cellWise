\name{cellMap}
\alias{cellMap}

\title{
Draw a cellmap
}

\description{
This function draws a cellmap, possibly of a subset of rows and columns of the data,
and possibly combining cells into blocks. A cellmap shows which cells are missing and which one are outlying, marking them in red for unusually large cell values and in blue for unusually low cell values. When cells are combined into blocks, the final color is the average of the colors in the individual cells. 
}

\usage{
cellMap(D, R, indcells, indrows, showVals = FALSE, xlabels = "",
        ylabels = "", mTitle = "",xtitle = "", ytitle = "",
        xshowindex = NULL, yshowindex = NULL, xblocksize = 1,
        yblocksize = 1, autolabel = TRUE, anglex = 90, sizexy = 1.1,
        hjustXlabels = 1,hjustYlabels = 1)
}

\arguments{
  \item{D}{
  The data matrix (required input argument).
  }
  \item{R}{
  Matrix of cell residuals (required input argument).
  }
  \item{indcells}{
  Indices of outlying cells (required input argument).
  }
  \item{indrows}{
  Indices of outlying rows (required input argument).
  }
  \item{showVals}{
  Whether to show the entries of D in the cellmap. Defaults to \code{FALSE}.
  }
  \item{xlabels}{
  Labels for the x-axis.
  }
  \item{ylabels}{
  Labels for the y-axis.
  }
  \item{mTitle}{
  Main title of the cellMap. 
  }
  \item{xtitle}{
  Title for the x-axis.
  }
  \item{ytitle}{
  Title for the y-axis.
  }
  \item{xshowindex}{
  Indices of the cells that will be shown, in the x direction.
  }
  \item{yshowindex}{
  Indices of the cells that will be shown, in the y direction.
  }
  \item{xblocksize}{
  Size of combination blocks in the \eqn{x} direction. Defaults to \eqn{1}.
  }
  \item{yblocksize}{
  Size of combination blocks in the \eqn{y} direction. Defaults to \eqn{1}.
  }
  \item{autolabel}{
  Automatically combines labels of cells in blocks.  
  If \code{FALSE}, you must provide the final xlabels and/or ylabels.
  Defaults to \code{TRUE}.
  }
  \item{anglex}{
  Angle of the labels on the x-axis. Defaults to \eqn{90}.
  }
  \item{sizexy}{
  Size of title for x-axis and y-axis. Defaults to \eqn{1.1}.
  }
  \item{hjustXlabels}{
  Adjust x-labels: 0=left, 0.5=centered, 1=right. Defaults to \eqn{1}.
  }
  \item{hjustYlabels}{
  Adjust y-labels: 0=left, 0.5=centered, 1=right. Defaults to \eqn{1}.
  }
}


\references{
 Rousseeuw P.J., Van den Bossche W. (2016). Detecting deviating data cells.
 \emph{arXiv:1601.07251}
}

\author{
 Rousseeuw P.J., Van den Bossche W. 
}

\seealso{
\code{\link{DetectDeviatingCells}}
}

\examples{
# For examples of the cellmap, we refer to the vignette:
vignette("DDC_examples")
}

