VonNeumann <- function(n, p=1, graine)
{
  x <-  rep(graine,n*p+1)
  for(i in 2:(n*p+1))
  {
    numbers <- strsplit(format(x[i-1]^2,scientific=FALSE),'')[[1]]
    while(length(numbers)>4){ 
        numbers <- numbers[2:(length(numbers)-1)] 
    }
    x[i] <- as.numeric(numbers)%*%(10^seq(length(numbers)-1,0,-1))
  }
  x <- matrix(x[2:(n*p+1)],nrow=n,ncol=p)
  return(x)
}


MersenneTwister <- function(n, p=1, graine)
{
  set.seed(graine,kind='Mersenne-Twister')
  x <- sample.int(2^32-1,n*p)
  x <- matrix(x,nrow=n,ncol=p)
  return(x)
}


binary <- function(x)
{
  if((x<2^31)&(x>=0))
    return( as.integer(intToBits(as.integer(x))) )
  else{
    if((x<2^32)&(x>0))
      return( c(binary(x-2^31)[1:31], 1) )
    else{
      cat('Erreur dans binary : le nombre etudie n est pas un entier positif en 32 bits.\n')
      return(c())
    }
  }
}

RANDU <- function(n=1,k = 10,graine)
{
  x <-  rep(graine,k*n+1)
  for (i in 2:(k*n+1)) {
    x[i] <- (65539*x[i-1])%%(2^31)
  }
  x <- matrix(x[2:(k*n+1)],nrow=n,ncol=k)
  
  return(x)
}

StandardMinimal <- function(n=1,k = 10,graine)
{
  x <-  rep(graine,k*n+1)
  
  for (i in 2:(k*n+1)) {
    x[i] <- (16807*x[i-1])%% (2^31-1)
  }
  
  x <- matrix(x[2:(k*n+1)],nrow=n,ncol=k)
  
  return(x)
}

Frequency <- function(x, nb)
{
  for(i in 1:nb){
    if(x[i]==0)
      x[i] <- -1
  }
  
  S <- sum(x)
  
  sObs <- abs(S)/sqrt(nb)
  
  p <- 2*(1-pnorm(sObs,0,1))
  
  return(p)
}

Runs <- function(x, nb)
{
  #Pre-test
  S <- 0
  for(i in 1:nb){
    if(x[i]==1)
      S <- S+1
  }
  pi <- S/nb
  
  if(abs(pi-1/2)>=(2/sqrt(nb))){
    return(0.0)
  }
  
  #Test
  VnObs <- 0
  for(i in 2:nb-1){
    if(x[i+1]!=x[i])
      VnObs <- VnObs + 1
  }
  VnObs <- VnObs + 1
  frac <- (abs(VnObs-2*nb*pi*(1-pi))/(2*sqrt(nb)*pi*(1-pi)))
  P <- 2*(1-pnorm(frac,0,1))
  
  return(P)
  
  
  
}