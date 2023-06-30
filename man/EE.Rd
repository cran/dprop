\name{Exponentiated exponential distribution}
\alias{d_EE}
\title{Compute the distributional properties of the exponentiated exponential distribution}
\usage{
d_EE(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the exponentiated exponential distribution.
}

\arguments{
\item{alpha}{
The strictly positive scale parameter of the exponential distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the  exponentiated exponential distribution (\eqn{\beta > 0}).
}
}

\details{The following is the probability density function of the  exponentiated exponential distribution:
\deqn{
f(x)=\alpha\beta e^{-\alpha x}\left(1-e^{-\alpha x}\right)^{\beta-1},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
d_EE gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the exponentiated exponential distribution.
}
\references{
Nadarajah, S. (2011). The exponentiated exponential distribution: a survey. AStA Advances in Statistical Analysis, 95, 219-251.

Gupta, R. D., & Kundu, D. (2007). Generalized exponential distribution: Existing results and some recent developments. Journal of Statistical Planning and Inference, 137(11), 3537-3547.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_EW}, \link{d_wei}, \link{d_exp}
}
}

\examples{
d_EE(5,2)

}
