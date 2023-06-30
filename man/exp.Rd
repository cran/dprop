\name{Exponential distribution}
\alias{d_exp}
\title{ Compute the distributional properties of the exponential distribution}
\usage{
d_exp(alpha)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the exponential distribution.
}

\arguments{
\item{alpha}{
The strictly positive scale parameter of the exponential distribution (\eqn{\alpha > 0}).
}
}

\details{The following is the probability density function of the exponential distribution:
\deqn{
f(x)=\alpha e^{-\alpha x},
}
where \eqn{x > 0} and \eqn{\alpha > 0}.
}
\value{
d_exp gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the exponential distribution.
}
\references{
Balakrishnan, K. (2019). Exponential distribution: theory, methods and applications. Routledge.

Singh, A. K. (1997). The exponential distribution-theory, methods and applications, Technometrics,
39(3), 341-341.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_wei}, \link{d_EE}
}
}

\examples{
d_exp(2)
}
