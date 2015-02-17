\name{getBC}
\alias{getBC}
\alias{getBC-method}
\title{
Get the details of the bicluster
}
\description{
\code{getBC} is a function to get the bicluster of \code{QUBIC}. \cr
\code{getBC} will return the gene names and condiction names of the bicluster you give. \cr \cr
thisBC = getBC(results, \ldots); \cr \cr
}
\arguments{
\item{results}{The result of QUBIC.}
\item{numBC}{Which bicluster you want to see.}
}
\value{
The R function, \code{getBC} returns an object of \code{list}: \cr 
\tabular{ll}{
genes \tab The gene names in this bicluster. \cr 
conds \tab The condiction names in this bicluster. \cr 
num_gene \tab The number of genes in this bicluster. \cr 
num_cond \tab The number of condictions in this bicluster. \cr 
}
}
\details{
results = qubic(data, file = "rQUBIC", q = 0.06, c = 0.95, f = 1, k = 2, r = 1, o = 100, d = "F"); \cr 
thisBC = getBC(results, numBC = 1);
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
\code{\link{QUBIC-package}} \code{\link{QUBIC-method}} \cr
}
\examples{

data(toy);
# load a discrete data named toy

results = qubic(toy, file = "toy", d = "T");
# run qubic to deal with the discrete data

thisBC = getBC(results, numBC = 3);
# get the gene and condiction names of the 3rd bicluster.

}
