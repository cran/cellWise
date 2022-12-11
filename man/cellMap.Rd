\name{cellMap}
\alias{cellMap}

\title{
Draw a cellmap
}

\description{
This function draws a cellmap, possibly of a subset of rows and columns of the data,
and possibly combining cells into blocks. A cellmap shows which cells are missing and which ones are outlying, marking them in red for unusually large cell values and in blue for unusually low cell values. When cells are combined into blocks, the final color is the average of the colors in the individual cells. 
}

\usage{
cellMap(R, indcells = NULL, indrows = NULL, standOD = NULL, 
          showVals = NULL, D=NULL, rowlabels = NULL,
          columnlabels = NULL, mTitle = "cell map",
          rowtitle = "cases", columntitle = "variables", 
          showrows = NULL, showcolumns = NULL,
          nrowsinblock = 1, ncolumnsinblock = 1, 
          autolabel = TRUE, columnangle = 90, sizemain = 2,
          sizetitles = 1.1, adjustrowlabels = 1, 
          adjustcolumnlabels = 1, colContrast = 1, outlyingGrad = TRUE, 
          darkestColor = sqrt(qchisq(0.999, 1)), drawCircles = TRUE) 
}

\arguments{
  \item{R}{
  Matrix of standardized residuals of the cells (required input argument). After running \code{\link{DDC}}, \code{\link{DDCpredict}},  \code{\link{MacroPCA}} or \code{\link{MacroPCApredict}} this is typically their value \code{$stdResid}.
  }
  \item{indcells}{
  Indices of outlying cells. Defaults to \code{NULL}, which indicates the cells for which
\eqn{|\code{R}| > \sqrt(qchisq(0.99,1))}.
  }
  \item{indrows}{
  Indices of outlying rows. If not \code{NULL}, the small circle to the right of the row is filled black if the row is in this list, and white otherwise. This gets overruled if \code{standOD} is not \code{NULL}.
  }
  \item{standOD}{
  Standardized Orthogonal Distance of each row. If not \code{NULL}, represents the outlyingness of each row by a shade of gray in the small circle to the right of the row. This color is white for \code{standOD} below \eqn{1}, and becomes fully black for \code{standOD} over \eqn{3}.
  }
  \item{showVals}{
  Takes the values \code{"D"}, \code{"R"} or \code{NULL} (the default). If \code{"R"} the numerical values of the residuals in \code{R} are shown in the cellmap. If \code{"D"}, the entries of the data matrix \code{D} are shown, provided the matrix \code{D} is being specified. If \code{NULL}, no entries are shown.
  }
  \item{D}{
  A matrix of data values, of the same dimensions as \code{R}. Default is \code{NULL}. \code{D} is only required when the data values are to be shown in the cellmap, by the option \code{showVals = "D"}. After running \code{\link{DDC}} or \code{\link{MacroPCA}}, \code{D} is typically their value \code{$remX}. After running \code{\link{DDCpredict}} or \code{\link{MacroPCApredict}} it is their argument \code{$newX}.
  }
  \item{rowlabels}{
  Labels of the rows. If \code{NULL}, row labels are taken from the argument \code{R}.
  }
  \item{columnlabels}{
  Labels of the columns. If \code{NULL}, column labels are taken from the argument \code{R}.
  }
  \item{mTitle}{
  Main title of the cellMap. Defaults to "cell map".
  }
  \item{sizemain}{
  Size of main title. Defaults to \eqn{2}.
  }
  \item{rowtitle}{
  Title for the rows. Defaults to "cases".
  }
  \item{columntitle}{
  Title for the columns. Defaults to "variables".
  }
  \item{sizetitles}{
  Size of row title and column title. Defaults to \eqn{1.1}.}
  \item{showrows}{
  Indices of the rows to be shown. Defaults to \code{NULL} which means all rows are shown.
  }
  \item{showcolumns}{
  Indices of the columns to be shown. Defaults to \code{NULL} which means all columns are shown.
  }
  \item{nrowsinblock}{
  How many rows are combined in a block. Defaults to \eqn{1}.
  }
  \item{ncolumnsinblock}{
  How many columns are combined in a block. Defaults to \eqn{1}.
  }
  \item{autolabel}{
  Automatically combines labels of cells in blocks.  
  If \code{FALSE}, you must provide the final columnlabels and/or rowlabels.
  Defaults to \code{TRUE}.
  }
  \item{columnangle}{
  Angle of the column labels. Defaults to \eqn{90} so the column labels are vertical.
  }
  \item{adjustrowlabels}{
  Adjust row labels: 0=left, 0.5=centered, 1=right. Defaults to \eqn{1}.
  }
  \item{adjustcolumnlabels}{
  Adjust column labels: 0=left, 0.5=centered, 1=right. Defaults to \eqn{1}.
  }
  \item{colContrast}{
 Parameter regulating the contrast of colors, should be in \eqn{[1,5]}. Defaults to \eqn{1}.
  }
  \item{outlyingGrad}{
  If \code{TRUE}, the color is gradually adjusted in function of the outlyingness. Defaults to \code{TRUE}.
  }
  \item{darkestColor}{
  Standardized residuals bigger than this will get the darkest color.
  }
  \item{drawCircles}{
  Whether or not to draw black circles indicating the outlying rows.
  }
}


\references{
Rousseeuw, P.J., Van den Bossche W. (2018). Detecting Deviating Data Cells. \emph{Technometrics}, \bold{60}(2), 135-145. \href{https://www.tandfonline.com/doi/full/10.1080/00401706.2017.1340909}{(link to open access pdf)}
}

\author{
 Rousseeuw P.J., Van den Bossche W. 
}

\seealso{
\code{\link{DDC}}
}

\examples{
# For examples of the cellmap, we refer to the vignette:
\dontrun{
vignette("DDC_examples")
}
}
