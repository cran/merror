alpha.beta.sigma <- function(x) {

  # Build alpha.beta.sigma matrix for use with cplot function from merror
  # Use with summary(fit)$parameters return from omx function
  
  k <- nrow(x)/3
  est <- x[,2]
  alpha <- c(0,est[(2*k+1):(3*k-1)])
  beta <- c(1,est[1:(k-1)])
  sigma <- sqrt(est[k:(2*k-1)])
  rbind(alpha,beta,sigma)
}
