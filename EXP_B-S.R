# ===================================== 
# German's Stochastic Volatility Model 
# Model: Explicit Black-Scholes 
# Date: April 27, 2013
# Author: Jacob Roth (jroth13@cmc.edu)
# =====================================

#   Asset:
S_0 <- 0
S_max <- 150
Ssteps <- M <- 100
S <- seq(S_0, S_max, length.out=Ssteps)
dS <- S[2] - S[1] # not S_max / Ssteps

#   Time:
t_0 <- 0
t_max <- capT <- 1
tsteps <- N <- 20000
t <- seq(t_0, t_max, length.out=tsteps)
tau = capT-t
dt <- t[2] - t[1]

#   Exogenous Parameters:
K = 75
sigma = .25
r_yr = .1
r = .1


#  Results Constructor
R <- array(NA,dim=c(M,N))
R[,N] = pmax(S-K,0)
R[1,] = 0
R[M,] = S_max - K*exp(-r*tau)

# Define Coefficients A, B, and C:
A <- function(i){
  A_i <- (1/2) * dt * (-r*i + sigma^2 * i^2)
  return(A_i)
}

B <- function(i){
  B_i <- 1 - dt * (sigma^2 * i^2 + r)
  return(B_i)
}

C <- function(i){
  C_i <- (1/2) * dt * (r*i + sigma^2 * i^2)
  return(C_i)
}

for (col in (N):2){
  for (row in (M-1):2){
    R[row,col-1] = A(row) * R[row-1,col] + B(row) * R[row,col] + C(row) * R[row+1,col]
  }
}  


    
shading <- rev(palette(rainbow(100)))
jet.colors <- colorRampPalette(shading[30:100])

nbcol <- 1000
color <- jet.colors(nbcol)
zfacet <- R[-1,-1] #+ R[-1,-(M+1)] + R[-(N+1),-1] + R[-(N+1),-1]
facetcol <- cut(zfacet,nbcol)
persp(R,col=color[facetcol],phi=10,theta=-50,shade=.00010,ltheta=90)
