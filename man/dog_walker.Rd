\name{dog_walker}
\alias{dog_walker}
\docType{data}
\title{
  The dog walker dataset
}
\description{
  A dataset containing the image sequence of a video. The sequence consists of 
  54 frames of 144 by 180 pixels pixels in Red/Geen/Blue (RGB) format.
}
\usage{data("dog_walker")}

\format{
  An array of dimensions \eqn{54 \times 144 \times 180 \times 3}.
}

\source{
 \url{http://www.wisdom.weizmann.ac.il/~vision/SpaceTimeActions.html}
}

\examples{
  data(dog_walker)
  # For more examples, we refer to the vignette:
vignette("Wrap_examples")
}
