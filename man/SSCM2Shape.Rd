\encoding{utf8}
\name{SSCM2Shape}
\alias{SSCM2Shape}

\title{
Calculation of the shape matrix
}
\description{
\code{SSCM2Shape} transforms the spatial sign covariance matrix of an elliptical distribution into its standardized shape matrix.
}
\usage{
SSCM2Shape(V,itermax=100,tol=10^(-10))
}
\arguments{
	\item{V}{ (required) p x p matrix representing the theoretical SSCM.}
	\item{tol}{ (optional) numeric, defines the stopping rule of the approximation procedure, see the help of \code{\link{evSSCM2evShape}} for details.}
	\item{itermax}{ (optional) numeric, defines the maximal number of iterations, see the help of \code{\link{evSSCM2evShape}} for details.}
}
\value{
  p x p symmetric numerical matrix, representing the shape matrix with trace 1, which corresponds to the spatial sign covariance matrix.
}
\details{
The calculation consists of three steps. First one calculates eigenvectors and eigenvalues of the SSCM matrix by the function \code{eigen}. Then one determines the eigenvalues of the related Shape matrix using the function \code{\link{evSSCM2evShape}}. Finally one expands the eigendecomposition consisting of the eigenvectors of the SSCM and the eigenvalues of the shape matrix. The resulting shape matrix is standardized to have a trace of 1. Note that this procedure only works for elliptical distributions.
}

\references{
\enc{Dürre}{Duerre}, A., Vogel, D., Fried, R. (2015): Spatial sign correlation, \emph{Journal of Multivariate Analyis}, vol. 135, 89--105. 
\href{http://arxiv.org/abs/1403.7635}{arvix 1403.7635}

\enc{Dürre}{Duerre}, A., Tyler, D. E., Vogel, D. (2016): On the eigenvalues of the spatial sign covariance matrix in more than two dimensions, to appear in: \emph{Statistics and Probability Letters}. \href{http://arxiv.org/abs/1512.02863}{arvix 1512.02863}

}
\seealso{
Calculating the theoretical shape from the theoretical SSCM \code{\link{SSCM2Shape}}

Calculating the eigenvalues of the SSCM \code{\link{evShape2evSSCM}}
}
\examples{

# defining a shape matrix with trace 1
V <- matrix(c(1,0.8,-0.2,0.8,1,0,-0.2,0,1),ncol=3)/3
V

# calculating the related SSCM
SSCM <- Shape2SSCM(V)

# recalculate the shape based on the SSCM
V2 <- SSCM2Shape(SSCM)
V2

# error is negligible
sum(abs(V-V2))
}
