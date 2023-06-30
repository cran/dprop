\name{Kumaraswamy exponential distribution}
\alias{d_kexp}
\title{Compute the distributional properties of the Kumaraswamy exponential distribution}
\usage{
d_kexp(lambda, a, b)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Kumaraswamy exponential distribution.

}

\arguments{
\item{a}{
The strictly positive shape parameter of the Kumaraswamy distribution (\eqn{a > 0}).
}
\item{b}{
The strictly positive shape parameter of the Kumaraswamy distribution (\eqn{b > 0}).
}
\item{lambda}{
The strictly positive parameter of the exponential distribution (\eqn{\lambda > 0}).
}


}

\details{The following is the probability density function of the Kumaraswamy exponential distribution:
\deqn{
f(x)=ab\lambda e^{-\lambda x}\left(1-e^{-\lambda x}\right)^{a-1}\left\{ 1-\left(1-e^{-\lambda x}\right)^{a}\right\} ^{b-1},
}
where \eqn{x > 0}, \eqn{a > 0}, \eqn{b > 0} and \eqn{\lambda > 0}.
}
\value{
d_kexp gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Kumaraswamy exponential distribution.


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
\code{\link{d_kburr}, \link{d_kum}
}
}

\examples{
d_kexp(0.2,1,1)

}
