#' @export
#' @import stats
d_wei <- function(alpha, beta) {
pdf_E <- function(x,alpha, beta) {
g=stats::dweibull(x,alpha, beta)
return(g)
}

E_wei <- function(alpha, beta, r) {
f <- function(alpha, beta, r, x) {
(x^(r))*(pdf_E(x,alpha, beta))
}
y=integrate(f,lower=0,upper=Inf,subdivisions=1000000,alpha=alpha,beta=beta, r=r)$value
return(y)
}
  m1 <- E_wei(alpha, beta, 1)
  m2 <- E_wei(alpha, beta, 2)
  m3 <- E_wei(alpha, beta, 3)
  m4 <- E_wei(alpha, beta, 4)
  mu1 <- m1
  mu2 <- m2-(m1)^2
  mu3 <- m3-3*(m1*m2)+2*(m1)^3
  mu4 <- m4-4*m3*m1+6*m2*(m1)^2-3*(m1)^4
  beta_1 <- (mu3^2)/(mu2^3)
  gamma_1 <- sqrt(beta_1)
  beta_2 <- mu4/(mu2)^2
  gamma_2 <- (beta_2)-3
  COV <- (sqrt(mu2)/mu1)*100



  med_wei <- qweibull(0.5,alpha, beta)

  qd_wei <- (qweibull(0.75,alpha, beta)-qweibull(0.25,alpha, beta))/2



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

  aux5 <- cbind(med_wei)
  colnames(aux5) = c("")

  aux6 <- cbind(qd_wei)
  colnames(aux6) = c("")


  list("First four ordinary moments"= round(aux,4), "First four central moments"= round(aux1,4), "mean and variance"=round(aux2,4), "Skewness and Kurtosis"=round(aux3,4), "Coefficient of variation"=round(aux4,4), "Median"=round(aux5,4),"Quartile deviation"=round(aux6,4))
}



