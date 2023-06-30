\name{Student's t distribution}
\alias{d_st}
\title{ Compute the distributional properties of the Student distribution}
\usage{
d_st(v)
}
\description{
Compute the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Student t distribution.
}

\arguments{
\item{v}{
The strictly positive parameter of the Student distribution (\eqn{v > 0}), it is also called a degree of freedom.
}
}

\details{The following is the probability density function of the Student t distribution:
\deqn{
f(x)=\frac{\Gamma(\frac{v+1}{2})}{\sqrt{v\pi}\Gamma(\frac{v}{2})}\left(1+\frac{x^{2}}{v}\right)^{-(v+1)/2},
}
where \eqn{x\in\left(-\infty,+\infty\right)} and \eqn{v > 0}.
}
\value{
d_st gives the first four ordinary moments, central moments, mean, and variance, Pearson's coefficient of skewness and kurtosis, coefficient of variation, median and quartile deviation based on the selected parametric values of the Student t distribution.
}
\references{
Yang, Z., Fang, K. T., & Kotz, S. (2007). On the Student's t-distribution and the t-statistic. Journal of Multivariate Analysis, 98(6), 1293-1304.

Ahsanullah, M., Kibria, B. G., & Shakil, M. (2014). Normal and Student's t distributions and their applications (Vol. 4). Paris, France: Atlantis Press.
}

\author{
Muhammad Imran.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{d_chi}
}
}

\examples{
d_st(6)
}
