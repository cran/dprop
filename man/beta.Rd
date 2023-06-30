\name{Beta distribution}
\alias{d_beta}
\title{Compute the distributional properties of the beta distribution}
\usage{
d_beta(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, variance, Pearson's coefficient of skewness, kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the beta distribution.
}

\arguments{
\item{alpha}{
The strictly positive shape parameter of the beta distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the beta distribution (\eqn{\beta > 0}).
}
}

\details{The following is the probability density function of the beta distribution:
\deqn{
f(x)=\frac{\Gamma\left(\alpha+\beta\right)}{\Gamma\left(\alpha\right)\Gamma\left(\beta\right)}x^{\alpha-1}\left(1-x\right)^{\beta-1},
}
where \eqn{0\leq x\leq1}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
d_beta gives the first four ordinary moments, central moments, mean, variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the beta distribution.
}
\references{
Gupta, A. K., & Nadarajah, S. (2004). Handbook of beta distribution and its applications. CRC Press.

Johnson, N. L., Kotz, S., & Balakrishnan, N. (1994). Beta distributions. Continuous univariate distributions. 2nd ed. New York, NY: John Wiley and Sons, 221-235.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_kum}
}
}

\examples{
d_beta(2,2)
}
