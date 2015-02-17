\name{Dtoy}
\alias{toy}
\docType{data}
\title{
A discrete data to test the argument -d
}
\description{
The toy is a matrix to test the discrete data. \cr
}
\usage{data(toy)}
\format{
  The format is:
 int [1:20, 1:12] 1 1 0 0 0 0 1 0 1 0 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:20] "G1" "G2" "G3" "G4" ...
  ..$ : chr [1:12] "C1" "C2" "C3" "C4" ...
}
\examples{
data(toy)
results = qubic(toy, d = "T");
}
\keyword{ datasets }
\keyword{ matrix }
\keyword{ discrete }
\seealso{
\code{\link{QUBIC-package}} \code{\link{QUBIC-method}} \cr
}