\name{Nakagami distribution}
\alias{d_naka}
\title{Compute the distributional properties of the Nakagami  distribution}
\usage{
d_naka(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Nakagami distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter of the Nakagami  distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter of the Nakagami  distribution (\eqn{\beta > 0}).
}

}

\details{The following is the probability density function of the Nakagami  distribution:
\deqn{
f(x)=\frac{2\alpha^{\alpha}}{\Gamma(\alpha)\beta^{\alpha}}x^{2\alpha-1}e^{-\frac{\alpha x^{2}}{\beta}},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
d_naka gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Nakagami distribution.
}
\references{
Schwartz, J., Godwin, R. T., & Giles, D. E. (2013). Improved maximum-likelihood estimation of the shape parameter in the Nakagami distribution. Journal of Statistical Computation and Simulation, 83(3), 434-445.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_gamma}
}
}

\examples{
d_naka(2,2)

}
