\name{Logistic distribution}
\alias{d_logis}
\title{Compute the distributional properties of the logistic distribution}
\usage{
d_logis(mu, sigma)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the logistic distribution.
}

\arguments{
\item{mu}{
Location parameter of the logistic distribution (\eqn{\mu\in\left(-\infty,+\infty\right)}).
}
\item{sigma}{
The strictly positive scale parameter of the logistic distribution (\eqn{\sigma > 0}).
}
}

\details{ The following is the probability density function of the logistic distribution:
\deqn{
f(x)=\frac{e^{-\frac{\left(x-\mu\right)}{\sigma}}}{\sigma\left(1+e^{-\frac{\left(x-\mu\right)}{\sigma}}\right)^{2}},
}
where \eqn{x\in\left(-\infty,+\infty\right)}, \eqn{\mu\in\left(-\infty,+\infty\right)} and \eqn{\sigma > 0}.
}
\value{
d_logis gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the logistic distribution.
}
\references{
Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). Continuous univariate distributions, Volume 2 (Vol. 289). John Wiley & Sons.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_lnormal}
}
}

\examples{
d_logis(4,0.2)
}
