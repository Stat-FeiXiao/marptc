gamma_0_with_var <-function(gS,var_alpha){
  kk <- length((gS))
  gamma_0 <- log(sum(gS))

  delta_g <- rep(1/sum(gS),kk)

  var_gamma <-  t(delta_g) %*% var_alpha %*% delta_g

  var_gamma<- diag(var_gamma)

  list(gamma_0=gamma_0, var_gamma= var_gamma)

}
