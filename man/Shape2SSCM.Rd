\encoding{utf8}
\name{Shape2SSCM}
\alias{Shape2SSCM}

\title{
Calculation of the Spatial Sign Covariance Matrix
}
\description{
\code{Shape2SSCM} transforms the theoretical shape matrix of an elliptical distribution into the spatial sign covariance matrix.
}
\usage{
Shape2SSCM(V)
}
\arguments{
	\item{V}{ (required) p x p matrix representing the theoretical shape matrix.}
}
\value{
  p x p symmetric numerical matrix, representing the spatial sign covariance matrix, which corresponds to the given shape matrix.
}
\details{
The calculation consists of three steps. First one calculates eigenvectors and eigenvalues of the shape matrix by the function \code{eigen}. Then one determines the related eigenvalues of the SSCM using the function \code{\link{evShape2evSSCM}} and finally one expands the resulting eigendecomposition consisting of the eigenvectors of the Shape matrix and the eigenvalues of the SSCM. Note that this procedure only works for elliptical distributions.
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
