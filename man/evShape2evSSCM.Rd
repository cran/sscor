\encoding{utf8}
\name{evShape2evSSCM}
\alias{evShape2evSSCM}

\title{
Calculation of the eigenvalues of the Spatial Sign Covariance Matrix
}
\description{
\code{evShape2evSSCM} transforms the eigenvalues of the shape matrix of an elliptical distribution into that of the spatial sign covariance matrix.
}
\usage{
evShape2evSSCM(evShape)
}
\arguments{
	\item{evShape}{ (required) p-dimensional numeric, representing the eigenvalues of the shape matrix.}
}
\value{
  p-dimensional numeric, representing the eigenvalues of the corresponding spatial sign covariance matrix.
}
\details{
The eigenvalues of the SSCM can be calculated from the eigenvalues of the shape matrix by numerical evaluation of onedimensional integrals, see Proposition 3 of \enc{Dürre}{Duerre}, Tyler, Vogel (2016). We use the substitution \deqn{x=\frac{1+t}{1-t}}{x=(1+t)/(1-t)} and Gaussian quadrature with Jacobi polynomials up to order 500 and \eqn{\beta=0}{beta=0} as well as \eqn{\alpha=p/2-1}{alpha=p/2-1}, see chapter 2.4 (iv) of Gautschi (1997) for details.

The nodes and weights of the Gauss-Jacobi-quadrature are originally computed by the \code{gaussquad} package and saved in the file \code{jacobiquad} for faster computation.
}
\references{
\enc{Dürre}{Duerre}, A., Vogel, D., Fried, R. (2015): Spatial sign correlation, \emph{Journal of Multivariate Analyis}, vol. 135, 89--105. 
\href{http://arxiv.org/abs/1403.7635}{arvix 1403.7635}

\enc{Dürre}{Duerre}, A., Tyler, D. E., Vogel, D. (2016): On the eigenvalues of the spatial sign covariance matrix in more than two dimensions, to appear in: \emph{Statistics and Probability Letters}. \href{http://arxiv.org/abs/1512.02863}{arvix 1512.02863}

Gautschi, W. (1997): \emph{Numerical Analysis - An Introduction}, \enc{Birkhäuser}{Birkhaeuser}, Basel.

Novomestky, F. (2013): gaussquad: Collection of functions for Gaussian quadrature. R package version 1.0-2.
}
\seealso{
Calculating the theoretical SSCM from the theoretical shape matrix \code{\link{Shape2SSCM}}
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
