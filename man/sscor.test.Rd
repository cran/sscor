\encoding{utf8}
\name{sscor.test}
\alias{sscor.test}

\title{
Correlation test based on spatial signs
}
\description{
Robust one-sample test and confidence interval for the correlation coefficient.
}

\usage{
sscor.test(x, y, rho0=0, alternative=c("two.sided","less","greater"),
conf.level=0.95, ...)
}

\arguments{
	\item{x,y}{(required) numeric vectors of observations, must have the same length.}
	\item{rho0}{(optional) correlation coefficient under the null hypothesis. The default is 0.}
	\item{alternative}{(optional) character string indicating the type of alternative to be tested. Must be one of \code{"two.sided"}, \code{"less"}, \code{"greater"}. The default is \code{"two-sided"}.}
	\item{conf.level}{(optional) confidence level. The default is 0.95.}
	\item{...}{optional arguments passed to sscor (such as location and scale estimates to be used).}
}	

\details{
The test is based on the spatial sign correlation (\enc{Dürre}{Duerre} et al. 2015), which is a highly robust correlation estimator, consistent for the generalized correlation coefficient under ellipticity. The confidence interval and the p-value are based on the asymptotic distribution after a variance-stabilizing transformation similar to Fisher's z-transform. They provide accurate approximations also for very small samples (\enc{Dürre}{Duerre} and Vogel, 2015).  The test is furthermore distribution-free within the elliptical model.  It has, e.g., the same power at the elliptical Cauchy distribution as at the multivariate Gaussian distribution.
}

\value{
A list with class \code{"htest"} containing the following values (similar to the output of \code{cor.test}):
	\item{statistic}{the value of the test statistic. Under the null, the test statistic is (asymptotically) standard normal.}
	\item{p.value}{the p-value of the test.}
	\item{estimate}{the estimated spatial sign correlation.}
	\item{null.value}{the true correlation under the null hypothesis.}
	\item{alternative}{a character string describing the alternative hypothesis.}
	\item{method}{a characters string indicating the choosen correlation estimator. Currently only the spatial sign correlation is implemented.}
	\item{data.name}{a character giving the names of the data.}
	\item{conf.int}{confidence interval for the correlation coefficient.}
}



\references{
\enc{Dürre}{Duerre}, A., Vogel, D., Fried, R. (2015): Spatial sign correlation, \emph{Journal of Multivariate Analyis}, vol. 135, 89--105. 
\href{http://arxiv.org/abs/1403.7635}{arvix 1403.7635}

\enc{Dürre}{Duerre}, A., Vogel, D. (2015): Asymptotics of the two-stage spatial sign correlation, preprint. \href{http://arxiv.org/abs/1506.02578}{arxiv 1506.02578}

}
\seealso{
Classical correlation testing: \code{\link{cor.test}}.

For more information on the spatial sign correlation: \code{\link{sscor}}.
}

\examples{

set.seed(5)
require(mvtnorm)

# create bivariate shape matrix with correlation 0.5
sigma <- matrix(c(1,0.5,0.5,1),ncol=2)

# under normality, both tests behave similarly
data <- rmvnorm(100,c(0,0),sigma)
x <- data[,1]
y <- data[,2]

sscor.test(x,y)
cor.test(x,y)

# sscor.test also works at a Cauchy distribution
data <- rmvt(100,diag(1,2), df=1)
x <- data[,1]
y <- data[,2]

sscor.test(x,y)
cor.test(x,y)
}
