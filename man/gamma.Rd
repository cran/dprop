\name{Gamma distribution}
\alias{d_gamma}
\title{Compute the distributional properties of the gamma distribution}
\usage{
d_gamma(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the gamma distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter of the gamma distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter of the gamma distribution (\eqn{\beta > 0}).
}
}

\details{The following is the probability density function of the gamma distribution:
\deqn{
f(x)=\frac{\beta^{\alpha}}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
d_gamma the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the gamma distribution.
}
\references{
Burgin, T. A. (1975). The gamma distribution and inventory control. Journal of the Operational Research Society, 26(3), 507-525.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_wei}, \link{d_naka}
}
}

\examples{
d_gamma(2,2)
}
