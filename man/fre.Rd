\name{Frechet distribution}
\alias{d_fre}
\title{Compute the distributional properties of the Frechet distribution}
\usage{
d_fre(alpha, beta, zeta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Frechet distribution.
}

\arguments{
\item{alpha}{
The parameter of the Frechet  distribution (\eqn{\alpha>0}).
}
\item{beta}{
The parameter of the Frechet  distribution (\eqn{\beta\in\left(-\infty,+\infty\right)}).
}
\item{zeta}{
The parameter of the Frechet  distribution (\eqn{\zeta>0}).
}
}

\details{The following is the probability density function of the Frechet distribution:

\deqn{
f(x)=\frac{\alpha}{\zeta}\left(\frac{x-\beta}{\zeta}\right)^{-1-\alpha}e^{-(\frac{x-\beta}{\zeta})^{-\alpha},}
}
where \eqn{x>\beta}, \eqn{\alpha>0}, \eqn{\zeta>0} and \eqn{\beta\in\left(-\infty,+\infty\right)}. The Frechet distribution is also known as inverse Weibull distribution and special case of the generalized extreme value distribution.
}
\value{
d_fre gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Frechet distribution.
}
\references{
Abbas, K., & Tang, Y. (2015). Analysis of Frechet distribution using reference priors. Communications in Statistics-Theory and Methods, 44(14), 2945-2956.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_wei}
}
}

\examples{
d_fre(5,1,0.5)

}
