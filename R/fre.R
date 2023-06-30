#' @export
#' @import VaRES
d_fre<-function(alpha, beta, zeta) {
  pdf_E<-function(x, alpha, beta, zeta){
    g=extraDistr::dfrechet(x, alpha, beta, zeta)
    return(g)
  }

  frechet_mom<-function(alpha, beta, zeta, r){
    f<-function(alpha, beta, zeta, r, x){
      (x^(r))*(pdf_E(x, alpha, beta, zeta))
    }
    y=integrate(f,lower=beta,upper=Inf,subdivisions=1000000,alpha=alpha,beta=beta,zeta=zeta, r=r)$value
    return(y)
  }
  m1<-frechet_mom(alpha, beta, zeta, 1)
  m2<-frechet_mom(alpha, beta, zeta, 2)
  m3<-frechet_mom(alpha, beta, zeta, 3)
  m4<-frechet_mom(alpha, beta, zeta, 4)
  mu1<-m1
  mu2<-m2-(m1)^2
  mu3<-m3-3*(m1*m2)+2*(m1)^3
  mu4<-m4-4*m3*m1+6*m2*(m1)^2-3*(m1)^4
  beta_1 <- (mu3^2)/(mu2^3)
  gamma_1 <- sqrt(beta_1)
  beta_2 <- mu4/(mu2)^2
  gamma_2 <- (beta_2)-3
  COV <- (sqrt(mu2)/mu1)*100



  med_frechet<-extraDistr::qfrechet(0.5, alpha, beta, zeta)

  qd_frechet<-(extraDistr::qfrechet(0.75, alpha, beta, zeta)-extraDistr::qfrechet(0.25, alpha, beta, zeta))/2



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

  aux5 <- cbind(med_frechet)
  colnames(aux5) = c("")

  aux6 <- cbind(qd_frechet)
  colnames(aux6) = c("")


  list("First four ordinary moments"= round(aux,4), "First four central moments"= round(aux1,4), "mean and variance"=round(aux2,4), "Skewness and Kurtosis"=round(aux3,4), "Coefficient of variation"=round(aux4,4), "Median"=round(aux5,4),"Quartile deviation"=round(aux6,4))
}

