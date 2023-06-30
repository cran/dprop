#' @export
#' @import stats
d_st <- function(v) {
pdf_E <- function(x, v) {
    g=stats::dt(x,v)
    return(g)
  }

  st_mom <- function(v,r) {
    f<-function(v,r,x) {
      (x^(r))*(pdf_E(x, v))
    }
    y=stats::integrate(f,lower=-Inf,upper=Inf,subdivisions=1000000,v=v, r=r)$value
    return(y)
  }
  m1 <- st_mom(v, 1)
  m2 <- st_mom(v, 2)
  m3 <- st_mom(v, 3)
  m4 <- st_mom(v, 4)
  mu1 <- m1
  mu2 <- m2-(m1)^2
  mu3 <- m3-3*(m1*m2)+2*(m1)^3
  mu4 <- m4-4*m3*m1+6*m2*(m1)^2-3*(m1)^4
  beta_1 <- (mu3^2)/(mu2^3)
  gamma_1 <- sqrt(beta_1)
  beta_2 <- mu4/(mu2)^2
  gamma_2 <- (beta_2)-3
  COV <- (sqrt(mu2)/mu1)*100

  med_st<-qt(0.5, v)

  qd_st<-(qt(0.75, v)-qt(0.25, v))/2

  aux <- cbind (m1, m2, m3, m4)
  colnames(aux) = c("Ist", "2nd",  "3rd",  "4th")

  aux1 <- cbind (mu1, mu2, mu3, mu4)
  colnames(aux1) = c("Ist", "2nd",  "3rd",  "4th")

  aux2 <- cbind (mu1, mu2)
  colnames(aux2) = c("mean", "variance")

  aux3 <- cbind (gamma_1, gamma_2)
  colnames(aux3) = c("Skewness", "Kurtosis")

  aux4 <- cbind (COV)
  colnames(aux4) = c("")

  aux5 <- cbind (med_st)
  colnames(aux5) = c("")

  aux6 <- cbind (qd_st)
  colnames(aux6) = c("")

  list("First four ordinary moments"= round(aux,4), "First four central moments"= round(aux1,4), "mean and variance"=round(aux2,4), "Skewness and Kurtosis"=round(aux3,4), "Coefficient of variation"=round(aux4,4), "Median"=round(aux5,4),"Quartile deviation"=round(aux6,4))
}




