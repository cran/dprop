\name{Chi-squared distribution}
\alias{d_chi}
\title{ Compute the distributional properties of the Chi-squared distribution}
\usage{
d_chi(n)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the (non-central) Chi-squared distribution.
}

\arguments{
\item{n}{
It is a degree of freedom and the positive parameter of the Chi-squared distribution (\eqn{n > 0}).
}
}

\details{The following is the probability density function of the (non-central) Chi-squared distribution:
\deqn{
f(x)=\frac{1}{2^{\frac{n}{2}}\Gamma(\frac{n}{2})}x^{\frac{n}{2}-1}e^{-\frac{x}{2}},
}
where \eqn{x > 0} and \eqn{n > 0}.
}
\value{
d_chi gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the (non-central) Chi-squared distribution.
}
\references{
Ding, C. G. (1992). Algorithm AS275: computing the non-central chi-squared distribution function. Journal of the Royal Statistical Society. Series C (Applied Statistics), 41(2), 478-482.

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
d_chi(2)
}
