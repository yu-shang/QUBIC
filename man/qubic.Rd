\name{QUBIC}
\alias{qubic}
\alias{QUBIC-method}
\alias{QUBIC}
\title{
Qualitative BiClustering Algorithm
}
\description{
\code{QUBIC} is a name short for QUalitative BIClustering algorithm, that can solve the biclustering problem in a more general form. \cr \cr 
The fuction \code{qubic} is the main function of the package. It calculates the bicluster in a data matrix, and it also discrete in a data matrix. \cr \cr
results = qubic(data, \ldots); \cr \cr
}
\arguments{
\item{data}{Data matrix.}
\item{file}{-file shoud be a string. The filename you want the QUBIC to write to. \cr default: "rQUBIC"}
\item{q}{-q should be (0, 0.5]. Use quantile discretization for continuous data.\cr default: 0.06 (see details in Method section in paper)}
\item{c}{-c should be (0.5, 1]. Consistency level of the block (0.5-1.0], the minimum ratio between the number of identical valid symbols in a column and the total number of rows in the output. \cr default: 0.95}
\item{f}{-f should be [0, 1.0]. Filtering overlapping blocks. \cr default: 1 (do not remove any blocks)}
\item{k}{-k should be >= 2. Minimum column width of the block. \cr default: 5\% of columns, minimum 2 columns}
\item{r}{-r should be a number. The number of ranks as which we treat the up(down)-regulated value when discretization. \cr default: 1}
\item{o}{-o should be a number. Number of blocks to report. \cr default: 100}
\item{d}{-d should be either ("True", "False"). Discrete data, where user should send their processed data to different value classes, see above. \cr default: False}
}
\value{
The R function, \code{qubic} returns an object of \code{list}: \cr 
\tabular{ll}{
Discrete \tab The discrete data of your given matrix. \cr 
ROWS_COLS \tab The rows and cols of the discrete data. \cr 
BC \tab The number of genes and conditions of every bicluster. \cr 
WriteToFile \tab The filename you have writen to. \cr 
ArgumentList \tab The parameters you have given to. \cr 
}
The R function, \code{qubic} writes servel files to your \code{Working Directory}, being the same as the C-code QUBIC does. (see details in paper) \cr

}
\details{
results = qubic(data, file = "rQUBIC", q = 0.06, c = 0.95, f = 1, k = 2, r = 1, o = 100, d = "F");
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
\code{\link{QUBIC-package}} \code{\link{getBC-method}} \code{\link{data}} \code{\link{toy}} \cr
}
\examples{

data(data);
# load the example matrix

results = qubic(data);
# run QUBIC with default parameters
# results is a variable of list with, 
#   the discrete data itself, 
#   the rows and cols of discrete data, 
#   the filename you want the QUBIC to write the results to, 
#   the argument list you sent to the QUBIC, 
# the more details of results, you will see in your Working Directory, 
# the name of that is what you sent to QUBIC (DEFAULT "rQUBIC"), 

thisBC = getBC(results, numBC = 1);
# get the gene and condiction names of the biggest, the 1st bicluster.

results = qubic(data, file = "AnotherExample");
# QUBIC will writing to your Working Directory with "AnotherExample" as file name

results = qubic(data, file = "Example", q = 0.06, c = 0.95, f = 1, k = 2, r = 1, o = 100, d = "F");
# run QUBIC with full parameters

data(toy);
# load a discrete data named toy

results = qubic(toy, file = "toy", d = "T");
# run qubic to deal with the discrete data
}
