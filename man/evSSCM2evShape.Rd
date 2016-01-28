\encoding{utf8}
\name{evSSCM2evShape}
\alias{evSSCM2evShape}

\title{
Calculation of the eigenvalues of the shape matrix
}
\description{
\code{evSSCM2evShape} transforms the eigenvalues of the SSCM of an elliptical distribution into that of the shape matrix.
}
\usage{
evSSCM2evShape(delta,tol=10^(-10),itermax=100)
}
\arguments{
	\item{delta}{ (required) p-dimensional numeric representing the eigenvalues of the SSCM.}
	\item{tol}{ (optional) numeric, defines the stopping rule of the approximation procedure, see details.}
	\item{itermax}{ (optional) numeric, defines the maximal number of iterations, see details.}
}
\value{
  p-dimensional numerical, representing the eigenvalues of the shape matrix. They are standardized to sum to 1.
}
\details{
The eigenvalues of the SSCM given that of the shape matrix can be calculated by evaluations of numerical integrals, see the help of \code{\link{evShape2evSSCM}} or \enc{Dürre}{Duerre}, Tyler, Vogel (2016). There is no expression for the inverse relationshop known. Though one can apply a fixed point iteration to get an approximation of the eigenvalues of the shape matrix. The iteration stops if either the maximal number of iterations is reached, which produces a warning, or if the error between the eigenvalues of the SSCM and the ones calculated from the actual fixed point iteration in L1 norm is smaller than the given tolerance. Since the mapping between the sets of eigenvalues is injective, see \enc{Dürre}{Duerre}, Tyler, Vogel (2016), this gives a reasonable approximation of the eigenvalues of the shape matrix.

}
\references{
\enc{Dürre}{Duerre}, A., Vogel, D., Fried, R. (2015): Spatial sign correlation, \emph{Journal of Multivariate Analyis}, vol. 135, 89--105. 
\href{http://arxiv.org/abs/1403.7635}{arvix 1403.7635}

\enc{Dürre}{Duerre}, A., Tyler, D. E., Vogel, D. (2016): On the eigenvalues of the spatial sign covariance matrix in more than two dimensions, to appear in: \emph{Statistics and Probability Letters}. \href{http://arxiv.org/abs/1512.02863}{arvix 1512.02863}

}
\seealso{
Calculating the theoretical shape from the theoretical SSCM \code{\link{SSCM2Shape}}

Calculating the eigenvalues of the SSCM from the eigenvalues of the shape matrix \code{\link{evShape2evSSCM}}
}
\examples{


# defining eigenvalues of the shape matrix
evShape <- seq(from=0,to=1,by=0.1)

# standardized to have sum 1
evShape <- evShape/sum(evShape)

# calculating the related eigenvalues of the SSCM
evSSCM <- evShape2evSSCM(evShape)

plot(evShape,evSSCM)

# recalculate the eigenvalues of the shape matrix
evShape2 <- evSSCM2evShape(evSSCM)

# error is negligible
sum(abs(evShape-evShape2))
}
