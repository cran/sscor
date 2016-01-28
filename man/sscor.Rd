\encoding{utf8}
\name{sscor}
\alias{sscor}

\title{
Spatial sign correlation
}
\description{
\code{sscor} computes a robust correlation matrix estimate based on spatial signs, as described in \enc{Dürre}{Duerre} et al. (2015).
}
\usage{
sscor(X, location=c("2dim-median","1dim-median","pdim-median","mean"),
scale=c("mad","Qn","sd"), standardized=TRUE, pdim=FALSE, ...)
}
\arguments{
	\item{X}{ (required) p x n data matrix, number of colums is the dimension p and the number of rows is the number of observations n.}
	\item{location}{ (optional) either a p-dimensional numeric vector specifying the location or a character string indicating the location estimator to be used. Possible values are \code{"2dim-median"},\code{"1dim-median"},\code{"pdim-median"},\code{"mean"}. The default is \code{"2dim-median"}. See details below.}
	\item{scale}{ (optional) either a p-dimensional numeric vector specifying the p marginal scales or a character string indicating the scale estimator to be used. Possible values are \code{"mad"},\code{"Qn"},\code{"sd"}. The default is  \code{"mad"}. See details below.}
	\item{standardized}{ (optional) logical; indicating whether the data should be standardized by marginal scale estimates prior to computing the spatial sign correlation. The default is \code{TRUE}.}
	\item{pdim}{ (optional) logical; indicating whether the correlation matrix consists of pairwise correlation estimates or is estimated at once by the p-dimensional spatial sign correlation, see details.}
	\item{...}{ (optional) arguments passed to \code{\link{evSSCM2evShape}} if \code{pdim=TRUE}.}
}
\value{
  p x p symmetric numerical matrix, the diagonal entries are 1, the off-diagonal entries are the pairwise spatial sign correlation estimates.
}
\details{
The spatial sign correlation is a highly robust estimator of the correlation matrix. It is consistent under elliptical distributions for the generalized correlation matrix (derived from the shape matrix instead of the correlation matrix, i.e., it is also defined when second moments are not finite).

There are two possibilities to calculate this matrix, one can either estimate all pairwise correlations by the two-dimensional spatial sign correlation or calculate the whole matrix at once by the p-dimensional spatial sign correlation. Both approaches have advantages and disadvantages. The first method should be more robust, especially if only some components of the observations are corrupted. Furthermore the consistency transformation is explicitly known only for the bivariate spatial sign correlation, whereas one has to apply an approximation procedure for the p-dimensional one. Additional argments can be passed to this algorithm using the \code{...} argument, see the help page of \code{\link{SSCM2Shape}} for details. On the other hand, the p-dimensional spatial sign correlation is more efficient under the normal distribution and always yields a positive semidefinite estimation.

The correlation estimator is computed in three steps: the data is standardized marginally, i.e., each variable is divided by a scale estimate. (This step is optional, but recommended, and hence the default.)
Then, if \code{pdim=FALSE}, for each pair of variables the 2x2 spatial sign covariance matrix (SSCM) is computed, and then from the SSCM a univariate correlation estimate given by the formulas (5) and (6) in \enc{Dürre}{Duerre} et al. (2015). These pairwise correlation estimates are the off-diagonal elements of the returned matrix estimate.
Otherwise, if {pdim=TRUE}, the pxp SSCM is computed, and then from the SSCM an estimator of the correlation matrix, which is done by the function \code{\link{SSCM2Shape}}, see there for details.


Scale estimation:


The scale estimates may either be computed outside the function \code{sscor} and passed on to \code{sscor} as a p-variate numeric vector, or they may be computed by \code{sscor}, using one of the following options: 

\code{"mad"}: applies \code{\link{mad}} from the standard package \code{stats}.  This is the default.

\code{"Qn"}: applies \code{\link{Qn}} from the package \code{robustbase}.

\code{"sd"}: applies the standard deviation \code{\link{sd}}. 

Standardizing the data is recommended (and is hence done by default), particularly so if the marginal scales largly differ. In this case, estimation without prior marginal standardization may become inefficient.


Location estimation:

The SSCM requires a multivariate location estimate. The location may be computed outside the function \code{sscor} and the result passed on to \code{sscor} as a p-variate numeric vector. Alternatively it may be computed by \code{sscor}, using one of the following options: 

\code{"2dim-median"}: two-dimensional spatial median, individually for every 2x2 SSCM. This is the default if \code{pdim=FALSE}.

\code{"1dim-median"}: the usual, one-dimensional median applied component-wise.

\code{"pdim-median"}: the p-dimensional spatial median for all variables. This is the default if \code{pdim=TRUE}.

\code{"mean"}: the p-dimensional mean. In light of robustness, it is not recommended to use the mean.


There is no handling of missing values.
}
\references{
\enc{Dürre}{Duerre}, A., Vogel, D., Fried, R. (2015): Spatial sign correlation, \emph{Journal of Multivariate Analyis}, vol. 135, 89--105. 
\href{http://arxiv.org/abs/1403.7635}{arvix 1403.7635}

\enc{Dürre}{Duerre}, A., Vogel, D. (2016): Asymptotics of the two-stage spatial sign correlation, \emph{Journal of Multivariate Analyis}, vol. 144, 54--67. \href{http://arxiv.org/abs/1506.02578}{arxiv 1506.02578}

}
\seealso{
Ordinary, non-robust correlation matrix: \code{\link{cor}}.

A number of other robust correlation estimators are provided by the package \code{rrcov}.

Testing for spatial sign correlation: \code{\link{sscor.test}}.
}
\examples{

set.seed(5)
X <- cbind(rnorm(25),rnorm(25))
# X is a 25x2 matrix

# cor() and sscor() behave similar under normality
sscor(X)
cor(X)

# but behave differently in the presence of outliers.
X[1,] <- c(10,10)
sscor(X)
cor(X) 

}
