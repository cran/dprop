\name{Chen distribution}
\alias{d_chen}
\title{Compute the distributional properties of the Chen distribution}
\usage{
d_chen(k, c)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Chen distribution.
}

\arguments{
\item{k}{
The strictly positive shape parameter of the Chen distribution (\eqn{k > 0}).
}
\item{c}{
The strictly positive scale parameter of the Chen distribution (\eqn{c > 0}).
}

}

\details{The following is the probability density function of the Chen distribution:
\deqn{
f(x)=ckx^{k-1}e^{x^{k}}e^{c-ce^{x^{k}}},
}
where \eqn{x > 0}, \eqn{c > 0} and \eqn{k > 0}.
}
\value{
d_chen gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Chen distribution.
}
\references{
Chen, Z. (2000). A new two-parameter lifetime distribution with bathtub shape or increasing failure rate function. Statistics & Probability Letters, 49(2), 155–161.
}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_wei}, \link{d_EE}, \link{d_EW}
}
}

\examples{
d_chen(0.2,0.2)

}
