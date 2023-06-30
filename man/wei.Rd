\name{Weibull distribution}
\alias{d_wei}
\title{Compute the distributional properties of the Weibull distribution}
\usage{
d_wei(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Weibull distribution.
}

\arguments{
\item{alpha}{
The strictly positive scale parameter of the Weibull distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the Weibull distribution (\eqn{\beta > 0}).
}
}

\details{The following is the probability density function of the Weibull distribution:
\deqn{
f(x)=\frac{\beta}{\alpha}\left(\frac{x}{\alpha}\right)^{\beta-1}e^{-(\frac{x}{\alpha})^{\beta}},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
d_wei gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Weibull distribution.
}
\references{
Hallinan Jr, Arthur J. (1993). A review of the Weibull distribution. Journal of Quality Technology, 25(2), 85-93.
}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_EE}
}
}

\examples{
d_wei(2,2)
}
