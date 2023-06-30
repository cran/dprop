\name{Burr XII distribution}
\alias{d_burr}
\title{Compute the distributional properties of the Burr XII distribution}
\usage{
d_burr(k, c)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Burr XII distribution.
}

\arguments{
\item{k}{
The strictly positive shape parameter of the Burr XII distribution (\eqn{k > 0}).
}
\item{c}{
The strictly positive shape parameter of the Burr XII distribution (\eqn{c > 0}).
}

}

\details{The following is the probability density function of the Burr XII distribution:
\deqn{
f(x)=kcx^{c-1}\left(1+x^{c}\right)^{-k-1},
}
where \eqn{x > 0}, \eqn{c > 0} and \eqn{k > 0}.
}
\value{
d_burr gives the first four ordinary moments, central moments, mean, variance, Pearson's coefficient of skewness, kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Burr XII distribution.
}
\references{
Rodriguez, R. N. (1977). A guide to the Burr type XII distributions. Biometrika, 64(1), 129-134.

Zimmer, W. J., Keats, J. B., & Wang, F. K. (1998). The Burr XII distribution in reliability analysis. Journal of Quality Technology, 30(4), 386-394.
}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_kburr}
}
}

\examples{
d_burr(2,10)

}
