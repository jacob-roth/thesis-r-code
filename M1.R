
~~~~~~~~~~~~
# =============================================
# German's Stochastic Volatility Model 
# Model: M1
# Date: April 27, 2013
# Author: Jacob Roth (jroth13@cmc.edu)
# Expected Run-Time: ~45 minutes for stability
# =============================================


# Discretize Variables
# -----------------------------------------------------------------------------------------
#   Asset:
S_0 <-   .01                                   
S_max <- 150                                   
Ssteps <- L <- 50                              
S <- l <- seq(S_0, S_max, length.out=Ssteps)   
dS = S_max / (Ssteps-1)

#   Volatility:
v_0 <- 0
v_max <- 1
vsteps <- W <- 10
v <- w <- seq(v_0, v_max, length.out=vsteps)
dv <- v_max / (vsteps-1)

#   Time:
t_0 <- 0
t_max <- capT <- 3
tsteps <- N <- 100000
t <- n <- seq(t_max, t_0, length.out=tsteps)
tau = capT- t
dt <- t_max / (tsteps-1)

#   Model Parameters
rho = .1
K = 75
r = .1
mu = .1
sigma = .25
# -----------------------------------------------------------------------------------------


# Initialize Boundary Conditions in Solution Arrays
# -----------------------------------------------------------------------------------------
#   Stock/Volatility Grid:
Uold <- array(NA,dim=c(L,W),dimnames=c("asset","vol"))
Unew <- array(NA,dim=c(L,W),dimnames=c("asset","vol"))
R <- array(NA,dim=c(L,W,N))


#   Payoff at Maturity (t = capT):
for (x in 1:L){
  for (y in 1:W){
    Uold[x,y] = pmax(S[x] - K, 0)
  }
}
#   t = capT BC
R[,,N] = Uold 

#   S = 0 BC
R[1,,] = 0                                                        

#   S --> inf BC
for (col in 1:W){
  for (row in 1:(N-1)){
    R[L,col,row] = pmax(S_max - K,0)
  } 
}

#   vol --> inf BC
for (col in 1:L){
  for (row in 1:(N-1)){
    R[col,W,row] = S[col]*exp(-r*t[row])
  } 
}                                                    

# Calculate Mesh
# -----------------------------------------------------------------------------------------
for (n_time in N:2){
  for (col in (W-1):2){
    for (row in (L-1):2){
#       # vol = 0 BC:
#       R[row,1,n_time-1] = dt/dv * mu*(R[row,1+1,n_time] 
#                         - R[row,1,n_time]) + 5*dt/dv^2 * sigma* (2*R[row,1,n_time] 
#                         - 5*R[row,1+1,n_time] + 4*R[row,1+2,n_time] 
#                         - R[row,1+3,n_time])
      # Interior Nodes:
      R[row,col,n_time-1] = dt*(mu*((R[row,col+1,n_time] - R[row,col-1,n_time])/dv) 
      + (1/2)*dv*col*row^2*((R[row+1,col,n_time]-2*R[row,col,n_time]+R[row-1,col,n_time])) 
      + (1/2)*sigma^2* ((R[row,col+1,n_time]-2*R[row,col,n_time]+R[row,col-1,n_time])/(dv^2)) 
      + rho*sigma*sqrt(dv)*sqrt(col)*row*((R[row+1,col+1,n_time]-R[row-1,col+1,n_time]-
                                  R[row+1,col-1,n_time]+R[row-1,col-1,n_time])/(4*dv))) 
      + R[row,col,n_time]
    }
  }
}
# =========================================================================================

# Display Resesults:
# -----------------------------------------------------------------------------------------
# REFERENCE: R[length,width,time]

#   Plotting Functions:
levelpersp <- function(x, y, z, colors=heat.colors, ...) {
  # getting the value of the midpoint
  zz <- (z[-1,-1] + z[-1,-ncol(z)] + z[-nrow(z),-1] + z[-nrow(z),-ncol(z)])/4
  # calculating the breaks
  breaks <- hist(zz, plot=FALSE)$breaks
  # cutting up zz
  cols <- rev(colors(length(breaks)-1))
  zzz <- cut(zz, breaks=breaks, labels=cols)
  # plotting
  persp(x, y, z, col=as.character(zzz), ...)
  # return breaks and colors for the legend
  list(breaks=breaks, colors=cols)
}

# Option Value at 20% Volatility
levelpersp(l,rev(t),R[,2,],theta=-30,phi=30,xlab="Asset Price",
           ylab="Time",zlab="Option Value",ticktype="detailed", 
           main="Option Value Over Time at 20% Volatility") 

# Option Value at t = 0.5
levelpersp(S[2:49],v[2:9],R[2:49,2:9,N/2],theta=-30,xlab="Asset Price",
           ylab="Time",zlab="Option Value",ticktype="detailed", 
           main="Option Value at t = 0.5")
