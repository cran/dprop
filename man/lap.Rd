\name{Laplace  distribution}
\alias{d_lap}
\title{Compute the distributional properties of the Laplace or double exponential distribution}
\usage{
d_lap(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Laplace distribution.
}

\arguments{
\item{alpha}{
Location parameter of the Laplace  distribution (\eqn{\alpha\in\left(-\infty,+\infty\right)}).
}
\item{beta}{
The strictly positive scale parameter of the Laplace  distribution (\eqn{\beta > 0}).
}
}

\details{The following is the probability density function of the Laplace  distribution:
\deqn{
f(x)=\frac{1}{2\beta}e^{\frac{-|x-\alpha|}{\beta}},
}
where \eqn{x\in\left(-\infty,+\infty\right)}, \eqn{\alpha\in\left(-\infty,+\infty\right)} and \eqn{\beta > 0}.
}
\value{
d_lap gives the first four ordinary moments, central moments, mean, variance, Pearson's coefficient of skewness, kurtosis, coefficient of variation, median and quartile deviation at some parametric values based on the Laplace  distribution.
}
\references{
Cordeiro, G. M., & Lemonte, A. J. (2011). The beta Laplace distribution. Statistics & Probability Letters, 81(8), 973-982.

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
d_lap(2,4)
}
