## test behaviour of ccr correction and slope of regression

## sample through SNR, nens, and nn space
SNRs <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)
nenses <- c(1,2,5,10,20,50,100,200,500)
nns <- c(5,10,20,30,50,100,200,500,1000, 2000, 5000, 1e5)
ntot <- length(SNRs)*length(nenses)*length(nns)

out <- data.frame(SNR=rep(SNRs, length(nenses)*length(nns)),
                 nens=rep(nenses, each=length(SNRs), length=ntot),
                 nn=rep(nns, each=length(SNRs)*length(nenses)))
out$rr <- out$ll <- out$corr <- NA
for (SNR in SNRs){
  for (nens in nenses){
    signal <- rnorm(max(nns))
    x <- rnorm(max(nns)) + sqrt(SNR) * signal
    y <- rnorm(max(nns), sd=sqrt(1 / nens)) + sqrt(SNR) * signal
    for (nn in nns){
      i <- which(out$SNR == SNR & out$nens == nens & out$nn == nn)
      sig_x <- sqrt(mean((x[1:nn] - mean(x[1:nn]))**2))
      sig_y <- sqrt(mean((y[1:nn] - mean(y[1:nn]))**2))
      out$rr[i] <- cor(x[1:nn], y[1:nn]) * sig_x / sig_y
      out$ll[i] <- coef(lm(x[1:nn] ~ y[1:nn]))[2]
      out$corr[i] <- cor(x[1:nn], y[1:nn])
      ## compute total least squares
      tsvd <- svd(cbind(y[1:nn], x[1:nn]))
      out$tls[i] <- -tsvd$v[1,2] / tsvd$v[2,2]
     }    
  }
}






