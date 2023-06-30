\name{Exponentiated Weibull distribution}
\alias{d_EW}
\title{Compute the distributional properties of the exponentiated Weibull distribution}
\usage{
d_EW(a, beta, zeta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the exponentiated Weibull distribution.
}

\arguments{
\item{a}{
The strictly positive shape parameter of the exponentiated Weibull distribution (\eqn{a > 0}).
}
\item{beta}{
The strictly positive scale parameter of the baseline Weibull distribution (\eqn{\beta > 0}).
}
\item{zeta}{
The strictly positive shape parameter of the baseline Weibull distribution (\eqn{\zeta > 0}).
}
}

\details{The following is the probability density function of the exponentiated Weibull distribution:
\deqn{
f(x)=a\zeta\beta^{-\zeta}x^{\zeta-1}e^{-\left(\frac{x}{\beta}\right)^{\zeta}}\left[1-e^{-\left(\frac{x}{\beta}\right)^{\zeta}}\right]^{a-1},
}
where \eqn{x > 0}, \eqn{a > 0}, \eqn{\beta > 0} and \eqn{\zeta > 0}.
}
\value{
d_EW gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the exponentiated Weibull distribution.
}
\references{
Nadarajah, S., Cordeiro, G. M., & Ortega, E. M. (2013). The exponentiated Weibull distribution: a survey. Statistical Papers, 54, 839-877.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_EE}, \link{d_wei}
}
}

\examples{
d_EW(1,1,0.5)

}
