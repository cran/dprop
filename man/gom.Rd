\name{Gompertz distribution}
\alias{d_gompertz}
\title{Compute the distributional properties of the Gompertz distribution}
\usage{
d_gompertz(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Gompertz distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter of the Gompertz  distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter of the Gompertz  distribution (\eqn{\beta > 0}).
}
}

\details{The following is the probability density function of the Gompertz distribution:
\deqn{
f(x)=\alpha e^{\beta x-\frac{\alpha}{\beta}\left(e^{\beta x}-1\right)},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
d_gompertz gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Gompertz distribution.
}
\references{
Soliman, A. A., Abd-Ellah, A. H., Abou-Elheggag, N. A., & Abd-Elmougod, G. A. (2012). Estimation of the parameters of life for Gompertz distribution using progressive first-failure censored data. Computational Statistics & Data Analysis, 56(8), 2471-2485.
}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_fre}
}
}

\examples{
d_gompertz(2,2)

}
