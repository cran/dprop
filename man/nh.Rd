\name{Exponential extension distribution}
\alias{d_nh}
\title{Compute the distributional properties of the exponential extension distribution}
\usage{
d_nh(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the exponential extension distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter of the exponential extension distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter of the exponential extension distribution (\eqn{\beta > 0}).
}
}

\details{The following is the probability density function of the exponential extension distribution:
\deqn{
f(x)=\alpha\beta(1+\alpha x)^{\beta-1}e^{1-(1+\alpha x)^{\beta}},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
d_nh gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the exponential extension distribution.
}
\references{
Nadarajah, S., & Haghighi, F. (2011). An extension of the exponential distribution. Statistics, 45(6), 543-558.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_exp}
}
}

\examples{
d_nh(0.5,1)

}
