#' @export
#' @import VaRES
d_kumnorm <- function(mu, sigma, a, b) {
  pdf_E <- function(x, mu, sigma, a, b) {
    g=VaRES::dkumnormal(x, mu, sigma, a, b)
    return(g)
  }

  knorm_mom <- function(mu, sigma, a, b, r) {
    f <- function(mu, sigma, a, b, r, x) {
      (x^(r))*(pdf_E(x, mu, sigma, a, b))
    }
    y=integrate(f,lower=-Inf,upper=Inf,subdivisions=1000000,mu=mu,sigma=sigma,a=a,b=b, r=r)$value
    return(y)
  }
  m1 <- knorm_mom(mu, sigma, a, b, 1)
  m2 <- knorm_mom(mu, sigma, a, b, 2)
  m3 <- knorm_mom(mu, sigma, a, b, 3)
  m4 <- knorm_mom(mu, sigma, a, b, 4)
  mu1 <- m1
  mu2 <- m2-(m1)^2
  mu3 <- m3-3*(m1*m2)+2*(m1)^3
  mu4 <- m4-4*m3*m1+6*m2*(m1)^2-3*(m1)^4
  alpha_1 <- (mu3^2)/(mu2^3)
  gamma_1 <- sqrt(alpha_1)
  alpha_2 <- mu4/(mu2)^2
  gamma_2 <- (alpha_2)-3
  COV <- (sqrt(mu2)/mu1)*100



  med_knorm<-VaRES::varkumnormal(0.5, mu, sigma, a, b)

  qd_knorm<-(VaRES::varkumnormal(0.75, mu, sigma, a, b)-VaRES::varkumnormal(0.25, mu, sigma, a, b))/2



  aux <- cbind(m1, m2, m3, m4)
  colnames(aux) = c("Ist", "2nd",  "3rd",  "4th")


  aux1 <- cbind(mu1, mu2, mu3, mu4)
  colnames(aux1) = c("Ist", "2nd",  "3rd",  "4th")

  aux2 <- cbind(mu1, mu2)
  colnames(aux2) = c("mean", "variance")

  aux3 <- cbind(gamma_1, gamma_2)
  colnames(aux3) = c("Skewness", "Kurtosis")

  aux4 <- cbind(COV)
  colnames(aux4) = c("")

  aux5 <- cbind(med_knorm)
  colnames(aux5) = c("")

  aux6 <- cbind(qd_knorm)
  colnames(aux6) = c("")


  list("First four ordinary moments"= round(aux,4), "First four central moments"= round(aux1,4), "mean and variance"=round(aux2,4), "Skewness and Kurtosis"=round(aux3,4), "Coefficient of variation"=round(aux4,4), "Median"=round(aux5,4),"Quartile deviation"=round(aux6,4))
}

