\name{Inverse-gamma distribution}
\alias{d_ingam}
\title{Compute the distributional properties of the inverse-gamma distribution}
\usage{
d_ingam(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the inverse-gamma distribution. }

\arguments{
\item{alpha}{
The strictly positive parameter of the inverse-gamma  distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter of the inverse-gamma  distribution (\eqn{\beta > 0}).
}
}

\details{The following is the probability density function of the inverse-gamma distribution:
\deqn{
f(x)=\frac{\beta^{\alpha}}{\Gamma(\alpha)}x^{-\alpha-1}e^{-\frac{\beta}{x}},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
d_ingam gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the inverse-gamma distribution.
}
\references{
Rivera, P. A., Calderín-Ojeda, E., Gallardo, D. I., & Gómez, H. W. (2021). A compound class of the inverse Gamma and power series distributions. Symmetry, 13(8), 1328.

Glen, A. G. (2017). On the inverse gamma as a survival distribution. Computational Probability Applications, 15-30.

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
d_ingam(5,2)

}
