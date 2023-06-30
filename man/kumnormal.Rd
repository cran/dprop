\name{Kumaraswamy normal  distribution}
\alias{d_kumnorm}
\title{Compute the distributional properties of the Kumaraswamy normal distribution}
\usage{
d_kumnorm(mu, sigma, a, b)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Kumaraswamy normal distribution.
}

\arguments{
\item{mu}{
The location parameter of the normal distribution (\eqn{\mu\in\left(-\infty,+\infty\right)).
}
}
\item{sigma}{
The strictly positive scale parameter of the normal distribution (\eqn{\sigma > 0}).
}
\item{a}{
The strictly positive shape parameter of the Kumaraswamy distribution (\eqn{a > 0}).
}
\item{b}{
The strictly positive shape parameter of the Kumaraswamy distribution (\eqn{b > 0}).
}
}


\details{The following is the probability density function of the Kumaraswamy normal distribution:
\deqn{
f(x)=\frac{ab}{\sigma}\phi\left(\frac{x-\mu}{\sigma}\right)\left[\Phi\left(\frac{x-\mu}{\sigma}\right)\right]^{a-1}\left[1-\Phi\left(\frac{x-\mu}{\sigma}\right)^{a}\right]^{b-1},
}
where \eqn{x\in\left(-\infty,+\infty\right)}, \eqn{\mu\in\left(-\infty,+\infty\right)}, \eqn{\sigma > 0}, \eqn{a > 0} and \eqn{b > 0}. The functions \eqn{\phi(.)} and \eqn{\Phi(.) }, denote the probability density function and cumulative distribution function of the standard normal variable, respectively.
}
\value{
d_kumnorm gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Kumaraswamy normal distribution.
}
\references{
Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions. Journal of Statistical Computation and Simulation, 81(7), 883-898.
}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_kburr}, \link{d_kexp}, \link{d_kum}
}
}

\examples{
d_kumnorm(0.2,0.2,2,2)

}
