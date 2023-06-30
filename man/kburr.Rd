\name{Kumaraswamy Burr XII distribution}
\alias{d_kburr}
\title{Compute the distributional properties of the Kumaraswamy Burr XII distribution}
\usage{
d_kburr(a, b, k, c)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Kumaraswamy Burr XII distribution.

}

\arguments{
\item{a}{
The strictly positive parameter of the Kumaraswamy distribution (\eqn{a > 0}).
}
\item{b}{
The strictly positive parameter of the Kumaraswamy distribution (\eqn{b > 0}).
}
\item{k}{
The strictly positive parameter of the Burr XII distribution (\eqn{k > 0}).
}
\item{c}{
The strictly positive parameter of the Burr XII distribution (\eqn{c > 0}).
}

}

\details{The following is the probability density function of the Kumaraswamy Burr XII distribution:
\deqn{
f(x)=\frac{abkcx^{c-1}}{\left(1+x^{c}\right)^{k+1}}\left[1-\left(1+x^{c}\right)^{-k}\right]^{a-1}\left\{ 1-\left[1-\left(1+x^{c}\right)^{-k}\right]^{a}\right\} ^{b-1},
}
where \eqn{x > 0}, \eqn{a > 0}, \eqn{b > 0}, \eqn{k > 0} and \eqn{c > 0}.
}
\value{
d_kburr gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Kumaraswamy Burr XII distribution.


}
\references{
Paranaiba, P. F., Ortega, E. M., Cordeiro, G. M., & Pascoa, M. A. D. (2013). The Kumaraswamy Burr XII distribution: theory and practice. Journal of Statistical Computation and Simulation, 83(11), 2117-2143.

}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_kum}, \link{d_kexp}
}
}

\examples{
d_kburr(1.5,1,1,7)

}
