\name{Lomax distribution}
\alias{d_lom}
\title{Compute the distributional properties of the Lomax distribution}
\usage{
d_lom(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Lomax distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter of the Lomax distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter of the Lomax distribution (\eqn{\beta > 0}).
}
}

\details{The following is the probability density function of the Lomax distribution:
\deqn{
f(x)=\frac{\alpha}{\beta}\left(1+\frac{x}{\beta}\right)^{-\alpha-1},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
d_lom gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Lomax distribution.
}
\references{
Abd-Elfattah, A. M., Alaboud, F. M., & Alharby, A. H. (2007). On sample size estimation for Lomax distribution. Australian Journal of Basic and Applied Sciences, 1(4), 373-378.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_gamma}
}
}

\examples{
d_lom(10,10)

}
