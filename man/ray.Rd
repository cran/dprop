\name{Rayleigh  distribution}
\alias{d_rayl}
\title{Compute the distributional properties of the Rayleigh distribution}
\usage{
d_rayl(alpha)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Rayleigh distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter of the Rayleigh distribution (\eqn{\alpha > 0}).
}

}

\details{The following is the probability density function of the Rayleigh distribution:
\deqn{
f(x)=\frac{x}{\alpha^{2}}e^{-\frac{x^{2}}{2\alpha^{2}}},
}
where \eqn{x > 0}, \eqn{\alpha > 0}.
}
\value{
d_rayl gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Rayleigh distribution.
}
\references{
Forbes, C., Evans, M. Hastings, N., & Peacock, B. (2011). Statistical Distributions. John Wiley &
Sons.
}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_wei}
}
}

\examples{
d_rayl(2)
}
