\name{Beta exponential distribution}
\alias{d_bexp}
\title{Compute the distributional properties of the beta exponential distribution}
\usage{
d_bexp(lambda, alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, variance, Pearson's coefficient of skewness, kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the beta exponential distribution.
}

\arguments{
\item{lambda}{
The strictly positive scale parameter of the exponential  distribution (\eqn{\lambda > 0}).
}
\item{alpha}{
The strictly positive shape parameter of the  beta distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the  beta  distribution (\eqn{\beta > 0}).
}
}

\details{The following is the probability density function of the  beta exponential distribution:
\deqn{
f(x)=\frac{\lambda e^{-\beta\lambda x}}{B(\alpha,\beta)}\left(1-e^{-\lambda x}\right)^{\alpha-1},
}
where \eqn{x > 0}, \eqn{\alpha > 0}, \eqn{\beta > 0} and \eqn{\lambda > 0}.
}
\value{
d_bexp gives the first four ordinary moments, central moments, mean, variance, Pearson's coefficient of skewness, kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the beta exponential distribution.
}
\references{
Nadarajah, S., & Kotz, S. (2006). The beta exponential distribution. Reliability Engineering & System Safety, 91(6), 689-697.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_beta}
}
}

\examples{
d_bexp(1,1,0.2)

}
