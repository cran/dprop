\name{Kumaraswamy distribution}
\alias{d_kum}
\title{Compute the distributional properties of the Kumaraswamy distribution}
\usage{
d_kum(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Kumaraswamy distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter of the Kumaraswamy distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter of the Kumaraswamy distribution (\eqn{\beta > 0}).
}
}

\details{ The following is the probability density function of the Kumaraswamy distribution:
\deqn{
f(x)=\alpha\beta x^{\alpha-1}\left(1-x^{a}\right)^{\beta-1},
}
where \eqn{0\leq x\leq1}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
d_kum gives the first four ordinary moments, central moments, mean, variance, Pearson's coefficient of skewness, kurtosis, coefficient of variation, median and quartile deviation at some parametric values based on the Kumaraswamy distribution.
}
\references{
El-Sherpieny, E. S. A., & Ahmed, M. A. (2014). On the kumaraswamy distribution. International Journal of Basic and Applied Sciences, 3(4), 372.

Mitnik, P. A. (2013). New properties of the Kumaraswamy distribution. Communications in Statistics-Theory and Methods, 42(5), 741-755.

Dey, S., Mazucheli, J., & Nadarajah, S. (2018). Kumaraswamy distribution: different methods of estimation. Computational and Applied Mathematics, 37, 2094-2111.


}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_beta}
}
}

\examples{
d_kum(2,2)
}
