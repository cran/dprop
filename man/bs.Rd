\name{Birnbaum-Saunders distribution}
\alias{d_bs}
\title{Compute the distributional properties of the Birnbaum-Saunders distribution}
\usage{
d_bs(v)
}
\description{
Compute the first four ordinary moments, central moments, mean, variance, Pearson's coefficient of skewness, kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Birnbaum-Saunders distribution.

}

\arguments{
\item{v}{
The strictly positive scale parameter of the Birnbaum-Saunders distribution (\eqn{v > 0}).
}

}

\details{The following is the probability density function of the Birnbaum-Saunders distribution:
\deqn{
f(x)=\frac{x^{0.5}+x^{-0.5}}{2vx}\phi\left(\frac{x^{0.5}-x^{-0.5}}{v}\right),
}
where \eqn{x > 0} and \eqn{v > 0}.
}
\value{
d_bs gives the first four ordinary moments, central moments, mean, variance, Pearson's coefficient of skewness, kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Birnbaum-Saunders distribution.
}
\references{
Chan, S., Nadarajah, S., & Afuecheta, E. (2016). An R package for value at risk and expected shortfall. Communications in Statistics Simulation and Computation, 45(9), 3416-3434.

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
d_bs(5)

}
