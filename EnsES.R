EnsES <- function(ens, obs, ...){
  nens <- ncol(ens)
  m <- rowMeans(ens)
  s <- apply(ens, 1, sd)
  g <- nens / (nens - 1) / (nens -2) * rowSums(((ens - m) / s)**3)
  e <- (m - obs)
  ES <- (s**2 - e**2 - e*s*g)**2
  return(ES)
}

EnsESss <- function (ens, ens.ref, obs){
  stopifnot(is.matrix(ens), is.matrix(ens.ref), is.vector(obs), 
            length(obs) == nrow(ens), length(obs) == nrow(ens.ref))
  xmask <- apply(!is.na(ens), 1, any) & !is.na(obs) & apply(!is.na(ens.ref), 
                                                            1, any)
  if (all(!xmask)) {
    xout <- NA
  }
  else {
    xout <- 1 - mean(EnsES(ens, obs))/mean(EnsES(ens.ref,obs))
  }
  return(xout)
}
