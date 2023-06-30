\name{F distribution}
\alias{d_F}
\title{Compute the distributional properties of the F distribution}
\usage{
d_F(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the F distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter of the F distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter of the F distribution (\eqn{\beta > 0}).
}

}

\details{The following is the probability density function of the F distribution:
\deqn{
f(x)=\frac{1}{B(\frac{\alpha}{2},\frac{\beta}{2})}\left(\frac{\alpha}{\beta}\right)^{\frac{\alpha}{2}}x^{\frac{\alpha}{2}-1}\left(1+\frac{\alpha}{\beta}x\right)^{-\left(\frac{\alpha+\beta}{2}\right)},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
d_F gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the F distribution.
}
\references{
Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). Continuous univariate distributions, volume 2 (Vol. 289). John Wiley & Sons.
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
d_F(2,10)

}
