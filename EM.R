# ===================================== 
# German's Stochastic Volatility Model 
# Model: Euler-Muruyama MC Simulation
# Date: April 27, 2013
# Author: Jacob Roth (jroth13@cmc.edu)
# =====================================


# Model Parameters (see reference [1])
eta <- 0
rho <- .1
K1 <- 75
K2 <- 125
K <- 75
dt <- .001

St <- 75
Vt <- .3

sigma <- .3
Vsum <- 0
a <- .1
b <- .1

EM <- function(n=100,type=1){
  W1 <- rnorm(n)
  W2 <- rnorm(1)
  
  pathS <- rep(NA,n+1)
  pathV <- rep(NA,n+1)
  
  pathS[1] <- St
  pathV[1] <- Vt
  
  for(time in 2:n){
    Vsum <- Vsum + Vt
    b <- Vsum/time
    mu <- a*(b-Vt)
    sigma <- sigma*sqrt(Vt)
    
    St <- St + eta*St*sqrt(Vt)*dt + St*sqrt(Vt)*W1[time]*sqrt(dt)
    repeat{
      W2 <- rho*W1[time] + sqrt(1-rho^2)*rnorm(1)
      Vt <- abs(mu*dt + sigma*W2*sqrt(dt))
      if(Vt > 0){
        break
      }
    }
    pathV[time] <- Vt
    pathS[time] <- St
  }
  payoffCall <- pmax(pathS[n] - K,0)
  payoffBCS <- pmin(pmax(pathS[n] - K1, 0), (K2 - K1))
  return(c(payoffCall,payoffBCS)[type])
}

EMprice <- function(trials){
  results <- replicate(trials,svpdeMC(100))
  mean <- mean(results)
  sd <- sd(results)/sqrt(trials)
  return(c(mean,sd))
}


