#' @export
#' @import VaRES
d_kburr <- function(a, b, k, c) {
  pdf_E <- function(x, a, b, k, c) {
    g=VaRES::dkumburr7(x, a, b, k, c)
    return(g)
  }

  kburr_mom <- function(a, b, k, c, r) {
    f <- function(a, b, k, c, r, x) {
      (x^(r))*(pdf_E(x, a, b, k, c))
    }
    y=integrate(f,lower=0,upper=Inf,subdivisions=1000000,a=a,b=b,k=k,c=c, r=r)$value
    return(y)
  }
  m1 <- kburr_mom(a, b, k, c, 1)
  m2 <- kburr_mom(a, b, k, c, 2)
  m3 <- kburr_mom(a, b, k, c, 3)
  m4 <- kburr_mom(a, b, k, c, 4)
  mu1 <- m1
  mu2 <- m2-(m1)^2
  mu3 <- m3-3*(m1*m2)+2*(m1)^3
  mu4 <- m4-4*m3*m1+6*m2*(m1)^2-3*(m1)^4
  alpha_1 <- (mu3^2)/(mu2^3)
  gamma_1 <- sqrt(alpha_1)
  alpha_2 <- mu4/(mu2)^2
  gamma_2 <- (alpha_2)-3
  COV <- (sqrt(mu2)/mu1)*100



  med_kburr <- VaRES::varkumburr7(0.5, a, b, k, c)

  qd_kburr <- (VaRES::varkumburr7(0.75, a, b, k, c)-VaRES::varkumburr7(0.25, a, b, k, c))/2



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

  aux5 <- cbind(med_kburr)
  colnames(aux5) = c("")

  aux6 <- cbind(qd_kburr)
  colnames(aux6) = c("")


  list("First four ordinary moments"= round(aux,4), "First four central moments"= round(aux1,4), "mean and variance"=round(aux2,4), "Skewness and Kurtosis"=round(aux3,4), "Coefficient of variation"=round(aux4,4), "Median"=round(aux5,4),"Quartile deviation"=round(aux6,4))
}


