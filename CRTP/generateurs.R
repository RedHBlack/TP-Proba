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
  cn <- binary(x[1])
  b <- cn[1:nb]
  for(i in 2:length(x)){
    cn <- binary(x[i])
    b <- c(b,cn[1:nb])
  }
  for(i in 1:length(b)){
    if(b[i]==0)
      b[i] <- -1
  }
  
  S <- sum(b)
  
  sObs <- abs(S)/sqrt(length(b))
  
  p <- 2*(1-pnorm(sObs,0,1))
  
  return(p)
}

Runs <- function(x, nb)
{
  cn <- binary(x[1])
  b <- cn[1:nb]
  for(i in 2:length(x)){
    cn <- binary(x[i])
    b <- c(b,cn[1:nb])
  }
  #Pre-test
  S <- sum(b)
  pi <- S/length(b)
  
  if(abs(pi-(1/2))>=(2/sqrt(length(b)))){
    return(0.0)
  }
  
  #Test
  VnObs <- 0
  for(i in 1:(length(b)-1)){
    if(b[i+1]!=b[i])
      VnObs <- VnObs + 1
  }
  VnObs <- VnObs + 1
  frac <- (abs(VnObs-2*length(b)*pi*(1-pi))/(2*sqrt(length(b))*pi*(1-pi)))
  P <- 2*(1-pnorm(frac,0,1))
  
  return(P)
}

FileMM1 <- function(lambda, mu, Duree){
  arrivee <- vector()
  depart <- vector()
  attente <- vector()
  T <- 0
  i <- 1
  while(T<=Duree){
    T <- T + rexp(1,lambda)
    if(T <= Duree){
      arrivee[i] <- T
      i <- i + 1
    }
  }
  depart[1] <- arrivee[1] + rexp(1, mu)
  attente[1] <- 0
  i <- 2
  while ( i<= length(arrivee) &&arrivee[i]<= Duree){
    depart[i] <- max(arrivee[i], depart[i-1])+ rexp(1, mu)
    attente[i] <- max(arrivee[i], depart[i-1]) - arrivee[i]
    i <- i+1
  }
  retour <- list(arrivee, depart, attente)
  return(retour)
}

FileMM2 <- function(lambda, mu, Duree){
  arrivee <- vector()
  depart <- vector()
  T <- 0
  i <- 1
  while(T<=Duree){
    T <- T + rexp(1,lambda)
    if(T <= Duree){
      arrivee[i] <- T
      i <- i + 1
    }
  }
  i <- 3
  depart[1] <- arrivee[1] + rexp(1, mu)
  depart[2] <- arrivee[2] + rexp(1, mu)
  while ( i<= length(arrivee) && arrivee[i]<= Duree){
    depart[i] <- max(arrivee[i], min(depart[i-1], depart[i-2]))+ rexp(1, mu)
    i <- i+1
  }
  retour <- list(arrivee, depart)
  return(retour)
}

nbClients <- function(arrivee, depart, enAttente){
  defaultW <- getOption("warn")
  options(warn = -1)
  i <- 1
  j <- 1
  k <- 1
  NBpersonnes <- vector()
  NBattente <- vector()
  NBpersonnes[k] <- 0
  NBattente[k] <- 0
  heure <- vector()
  heure[k] <- 0
  k <- k+1
  while(i <= length(arrivee) || j <= length(depart) ){
    if(i <= length(arrivee) && arrivee[i] < min(depart, na.rm = TRUE)){
      NBpersonnes[k] <- NBpersonnes[k-1]+1
      
      heure[k] <- arrivee[i]
      i <- i + 1
    }else{
      NBpersonnes[k] <- NBpersonnes[k-1]-1
      heure[k] <- min(depart, na.rm = TRUE)
      depart <- depart[depart != heure[k]]
      j <- j + 1
    }
    if(NBpersonnes[k] != 0){
      NBattente[k] = NBpersonnes[k] -1
    }else{
      NBattente[k] = 0
    }
    k <- k +1
  }
  options(warn = defaultW)
  systeme <- list(heure, NBpersonnes, NBattente)
  return(systeme)
}
# Moyenne nb clients et moyenne temps resté
Moyennes <- function(nbClients, fileDepArr, Duree){
  moyNbClients <- vector()
  for(i in 1:length(nbClients[[2]])){
    if(i!=length(nbClients[[2]])){
      moyNbClients[i] <- nbClients[[2]][i]*(nbClients[[1]][i+1]-nbClients[[1]][i])
    }else{
      moyNbClients[i] <- nbClients[[2]][i]*(Duree-nbClients[[1]][i])
    }
  }
  moyNbClients <- sum(moyNbClients)/Duree
  tempsMoyen <- sum(fileDepArr[[2]]-fileDepArr[[1]])/length(fileDepArr[[1]])
  
  Na <- vector()
  for(i in 1:length(nbClients[[3]])){
    if(i!=length(nbClients[[3]])){
      Na[i] <- nbClients[[3]][i]*(nbClients[[1]][i+1]-nbClients[[1]][i])
    }else{
      Na[i] <- nbClients[[3]][i]*(Duree-nbClients[[1]][i])
    }
  }
  Na <- sum(Na)/Duree
  Wa <- mean(fileDepArr[[3]])
  
  moyennes <- list(moyNbClients, tempsMoyen, Na, Wa)
  return (moyennes)
}























































































