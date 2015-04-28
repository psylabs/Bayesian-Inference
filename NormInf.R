NormInf <- function(data, mu0, s20, k0, nu0){
  
  ##   GIVEN CONJUGATE NORM PRIOR & NORM SAMPLING MODEL         ##
  ##   COMPUTES POSTERIOR MEAN & CI AND PLOTS PRIOR/POSTERIOR   ##
  ##   FOR BOTH MU AND S2   ##
 
  n    <- length(data)
  s2   <- var(data)
  ybar <- mean(data)
  kn   <- k0 + n
  nun  <- nu0 + n
  
  mun <- k0 * mu0/kn + n * ybar/kn
  s2n <- (nu0*s20 + (n-1)*s2 + k0*n * (ybar - mu0)^2/kn) / nun
  
  # MC POSTERIOR DISTRIBUTIONS #
  
  s2.mc <- 1/rgamma(10000, nun/2, nun*s2n/2)
  mu.mc <- rnorm(10000, mun, sqrt(s2.mc/kn))        #Why do we divide s2/kn
  q <- data.frame(mu = quantile(mu.mc,c(.025,.5,.975)), s2 = quantile(s2.mc, c(.025,.5,.975))) 
  row.names(q) <- c("2.5%", "50%", "97.5%")
  print(q)
  # PLOT 4 DISTRIBUTIONS PRIOR & POST FOR MU AND S2
  
  par(mfrow = c(2,1))
  plot(density(mu.mc), xlim = c(-2*min(mu.mc), +2*max(mu.mc)), xlab="mu", ylab="density", type="l",
       lwd=2, col=2, lty=2, bty="n", main ="mu")
  lines(density(rnorm(10000,mu0,sqrt(s20))), lwd=2, col=1, lty=1) 
  legend("topright", c("prior", "post"),
         lwd=2, col=1:2, lty=1:2, bty="n")    
  
  plot(density(s2.mc), xlim = c(-2*min(s2.mc), +2*max(s2.mc)), xlab="s2", ylab="density", type="l",
       lwd=2, col=2, lty=2, bty="n", main ="s2")
  s2.prior <- 1/rgamma(10000, nu0/2, nu0*s20/2)
  lines(density(s2.prior), lwd=2, col=1, lty=1)  #NOT SURE IF THIS IS A SENSIBLE PRIOR
  legend("topright", c("prior", "post"),
         lwd=2, col=1:2, lty=1:2, bty="n")    
  
  return (mu.mc)  
}

m1 <- NormInf(data1 <- scan("school1.txt"), mu0=5, s20=4, k0=1, nu0=2)
m2 <- NormInf(data1 <- scan("school2.txt"), mu0=5, s20=4, k0=1, nu0=2)
m3 <- NormInf(data1 <- scan("school3.txt"), mu0=5, s20=4, k0=1, nu0=2)

mx <-cbind(m1,m2,m3)
m <- as.list(data.frame(t(cbind(m1,m2,m3))))

# b
library(gtools)
x <- c(1,2,3)
perms <- permutations(n=3,r=3,v=x)

comp <- function(m,i,j,k){
  if (m[i]<m[j] && m[j]<m[k]) {return (TRUE)}
  else { return(FALSE)}
}

tot = 0
for (perm in split(perms,row(perms))){
  i = perm[1];j = perm[2];k = perm[3]
  message(sprintf("%s < %s < %s",i,j,k))
  a <-apply(mx,1,comp,i=i,j=j,k=k)
  print(mean(a))
  tot = tot + mean(a)
}


# why is mu50, not mun