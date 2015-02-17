\name{QUBIC-package}
\alias{QUBIC-package}
\docType{package}
\title{
QUBIC is a Qualitative BiClustering Algorithm for Expression Data Analysis
}
\description{
  Biclustering extends the traditional clustering techniques by attempting to find (all) subgroups of genes with similar expression patterns under to-be-identified subsets of experimental conditions when applied to gene expression data. Still the real power of this clustering strategy is yet to be fully realized due to the lack of effective and efficient algorithms for reliably solving the general biclustering problem. We report a QUalitative BIClustering algorithm (QUBIC) that can solve the biclustering problem in a more general form, compared to existing algorithms, through employing a combination of qualitative (or semi-quantitative) measures of gene expression data and a combinatorial optimization technique. One key unique feature of the QUBIC algorithm is that it can identify all statistically significant biclusters including biclusters with the so-called 'scaling patterns', a problem considered to be rather challenging; another key unique feature is that the algorithm solves such general biclustering problems very efficiently, capable of solving biclustering problems with tens of thousands of genes under up to thousands of conditions in a few minutes of the CPU time on a desktop computer. We have demonstrated a considerably improved biclustering performance by our algorithm compared to the existing algorithms on various benchmark sets and data sets of our own. 
 \cr \cr
Certain parts of the code uses open-source data structure library codes, including: \cr
- fib <\url{http://resnet.uoregon.edu/~gurney_j/jmpc/fib.html}>, copyright information in fib.c \cr
- Mark A. Weiss's data structure codes <\url{http://www.cs.fiu.edu/~weiss/}> \cr
 \cr
The main function of the package is qubic(), and the more details you can find in \code{\link{QUBIC-method}} \cr
}
\details{
\tabular{ll}{
Package: \tab QUBIC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-11\cr
Depends: \tab R (>= 3.0.3), Rcpp (>= 0.11.3)\cr
LinkingTo: \tab Rcpp\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
LazyData: \tab true\cr
}
}
\author{
Qin Ma (UGA) <\email{maqin@uga.edu}> \cr 
Yu Shang (JLU & UGA)<\email{yushang@uga.edu}> \cr 
Ying Xu (JLU & UGA) <\email{xyn@bmb.uga.edu}> \cr \cr 
Maintainer: Yu Shang (JLU & UGA)<\email{yushang@uga.edu}> \cr 
}
\references{
Li, Ma, et al. (2009) \emph{QUBIC: a qualitative biclustering algorithm for analyses of gene expression data} Nucleic Acids Research 07/2009 37(15):e101 \cr
\url{https://code.google.com/p/qubic-biclustering/} \cr\cr 
Zhou, Ma, et al. (2012) \emph{QServer: a biclustering server for prediction and assessment of co-expressed gene clusters} PLoS ONE 01/2012 7(3):e32660 \cr
}
\keyword{ cluster }
\keyword{ bicluster }
\seealso{
\code{\link{QUBIC-method}} \code{\link{getBC-method}} \cr
}
