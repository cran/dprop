\name{Gumbel distribution}
\alias{d_gumbel}
\title{Compute the distributional properties of the Gumbel distribution}
\usage{
d_gumbel(alpha, beta)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Gumbel distribution.
}
\arguments{
\item{alpha}{
Location parameter of the Gumbel distribution (\eqn{\alpha\in\left(-\infty,+\infty\right)}).
}
\item{beta}{
The strictly positive scale parameter of the Gumbel distribution (\eqn{\beta > 0}).
}
}
\details{The following is the probability density function of the Gumbel distribution:
\deqn{f(x)=\frac{1}{\beta}e^{-(z+e^{-z})},}
where \eqn{z=\frac{x-\alpha}{\beta}}, \eqn{x\in\left(-\infty,+\infty\right)}, \eqn{\alpha\in\left(-\infty,+\infty\right)} and \eqn{\beta > 0}.
}
\value{
d_gumbel gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Gumbel distribution.
}
\references{

Gomez, Y. M., Bolfarine, H., & Gomez, H. W. (2019). Gumbel distribution with heavy tails and applications to environmental data. Mathematics and Computers in Simulation, 157, 115-129.
}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_gompertz}, \link{d_fre}
}
}

\examples{
d_gumbel(1,2)
}
