\name{Normal distribution}
\alias{d_normal}
\title{Compute the distributional properties of the normal distribution}
\usage{
d_normal(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the normal distribution.

}

\arguments{
\item{alpha}{
Location parameter of the normal distribution (\eqn{\alpha\in\left(-\infty,+\infty\right)}).
}
\item{beta}{
The strictly positive scale parameter of the normal distribution (\eqn{\beta > 0}).
}
}

\details{The following is the probability density function of the normal distribution:
\deqn{
f(x)=\frac{1}{\beta\sqrt{2\pi}}e^{-0.5\left(\frac{x-\alpha}{\beta}\right)^{2}},
}
where \eqn{x\in\left(-\infty,+\infty\right)}, \eqn{\alpha\in\left(-\infty,+\infty\right)} and \eqn{\beta > 0}. The parameters \eqn{\alpha} and \eqn{\beta} represent the mean and standard deviation, respectively.
}
\value{
d_normal gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the normal distribution.


}
\references{
Patel, J. K., & Read, C. B. (1996). Handbook of the normal distribution (Vol. 150). CRC Press.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_lnormal}
}
}

\examples{
d_normal(4,0.2)
}
