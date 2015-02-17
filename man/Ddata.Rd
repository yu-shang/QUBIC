\name{Ddata}
\alias{data}
\docType{data}
\title{
A continuous data of Ecoli
}
\description{
The data is a matrix to test the continuous data. \cr
}
\usage{data(data)}
\format{
  The format is:
 num [1:1500, 1:466] 9.08 7.12 7.18 7.9 6.8 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:1500] "b4634" "b3241" "b3240" "b3243" ...
  ..$ : chr [1:466] "dinI_U_N0025" "dinP_U_N0025" "lexA_U_N0025" "lon_U_N0025" ...
}
\examples{
data(data)
results = qubic(data);
}
\keyword{ datasets }
\keyword{ matrix }
\keyword{ continuous }
\seealso{
\code{\link{QUBIC-package}} \code{\link{QUBIC-method}} \cr
}