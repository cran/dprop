\name{Log-normal distribution}
\alias{d_lnormal}
\title{Compute the distributional properties of the log-normal distribution}
\usage{
d_lnormal(mu, sigma)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the log-normal distribution.
}

\arguments{

\item{mu}{
The location parameter (\eqn{\mu\in\left(-\infty,+\infty\right)}).
}
\item{sigma}{
The strictly positive scale parameter of the log-normal distribution (\eqn{\sigma > 0}).
}
}

\details{The following is the probability density function of the log-normal distribution:
\deqn{
f(x)=\frac{1}{x\sigma\sqrt{2\pi}}e^{-\frac{\left(\log(x)-\mu\right)^{2}}{2\sigma^{2}}},
}
where \eqn{x > 0}, \eqn{\mu\in\left(-\infty,+\infty\right)} and \eqn{\sigma > 0}.
}
\value{
d_lnormal gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the log-normal distribution.
}
\references{
Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). Continuous Univariate Distributions, Volume 1, Chapter 14. Wiley, New York.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_normal}
}
}

\examples{
d_lnormal(1,0.5)
}
