#' @export
#' @import VaRES
d_EW <- function(a, beta, zeta) {
pdf_E <- function(x,a, beta, zeta) {
    g=VaRES::dexpweibull(x, a, beta, zeta)
    return(g)
  }

  EW_mom <- function(a, beta, zeta, r) {
    f <- function(a, beta, zeta, r, x) {
      (x^(r))*(pdf_E(x, a, beta, zeta))
    }
    y=integrate(f,lower=0,upper=Inf,subdivisions=1000000,a=a,beta=beta,zeta=zeta, r=r)$value
    return(y)
  }
  m1 <- EW_mom(a, beta, zeta, 1)
  m2 <- EW_mom(a, beta, zeta, 2)
  m3 <- EW_mom(a, beta, zeta, 3)
  m4 <- EW_mom(a, beta, zeta, 4)
  mu1 <- m1
  mu2 <- m2-(m1)^2
  mu3 <- m3-3*(m1*m2)+2*(m1)^3
  mu4 <- m4-4*m3*m1+6*m2*(m1)^2-3*(m1)^4
  beta_1 <- (mu3^2)/(mu2^3)
  gamma_1 <- sqrt(beta_1)
  beta_2 <- mu4/(mu2)^2
  gamma_2 <- (beta_2)-3
  COV <- (sqrt(mu2)/mu1)*100



  med_EW <- VaRES::varexpweibull(0.5, a, beta, zeta)

  qd_EW <- (VaRES::varexpweibull(0.75, a, beta, zeta)-VaRES::varexpweibull(0.25 , a, beta, zeta))/2



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

  aux5 <- cbind(med_EW)
  colnames(aux5) = c("")

  aux6 <- cbind(qd_EW)
  colnames(aux6) = c("")


  list("First four ordinary moments"= round(aux,4), "First four central moments"= round(aux1,4), "mean and variance"=round(aux2,4), "Skewness and Kurtosis"=round(aux3,4), "Coefficient of variation"=round(aux4,4), "Median"=round(aux5,4),"Quartile deviation"=round(aux6,4))
}


