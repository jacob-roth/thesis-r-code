~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ===================================== 
# German's Stochastic Volatility Model 
# Model: M2
# Date: April 27, 2013
# Author: Jacob Roth (jroth13@cmc.edu)
# =====================================


# Discretize Variables
# -----------------------------------
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
t_max <- capT <- 1
tsteps <- N <- 10000
t <- n <- seq(t_max, t_0, length.out=tsteps)
tau = capT- t
dt <- t_max / (tsteps-1)

#   #   (constant) Drift / Vol Functions:
#   mu <- function(t,v){
#     return(.1)
#   }
#   sigma <- function(t,v){
#     return(.25)
#   }
#   
rho = .1
K1 = 75
K2 = 125
r = 0
#mu = .1
#sigma = .25


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
    Uold[x,y] = pmin(pmax(S[x] - K1, 0), K2 - K1)
  }
}
#   Boundary Conditions

#   s = 0 BC
R[1,,] = 0                                                        

#   s --> inf BC
for (col in 1:W){
  for (row in 1:(N-1)){
    R[L,col,row] = pmin(pmax(S_max - K1*exp(-r*t[row]), 0), (K2-K1)*exp(-r*t[row]))
  }
}

#VERSION 2 --> wiht vol zero
for (col in 1:L){
  for (row in 1:(N-1)){
    if (S[col] >= K2){
      R[col,1,row] = (K2-K1)*exp(-r*t[row]) # t acts like ttm
    }
    if (S[col] < K2){
      R[col,1,row] = pmin(pmax(S[col]-K1,0),(K2-K1))
    }
    }
  }

#   vol --> inf BC
for (col in 1:L){
  for (row in 1:(N-1)){
    if (S[col] >= K2-K1){
      R[col,W,row] = (K2-K1)*exp(-r*t[row]) # t acts like ttm
    }
    if (S[col] < K2-K1){
      R[col,W,row] = S[col]*exp(-r*t[row])
    
  } 
}  
}

R[,,N] = Uold                                                   



# ---------------------------------------

Vsum <- 0
a <- 0.1

# Mesh Calculation 
# ---------------------------------------
for (n_time in N:2){
  for (col in (W-1):2){
    Vsum <- Vsum + col*dv
    b <- Vsum/n_time
    mu <- a*(b-dv*col)
    sigma <- sigma * sqrt(dv*col)
    for (row in (L-1):2){
      R[row,col,n_time-1] = dt*(mu*((R[row,col+1,n_time] - R[row,col-1,n_time])/dv) + 
                            (1/2)*dv*col*row^2* ((R[row+1,col,n_time]-2*R[row,col,n_time]+R[row-1,col,n_time])) + 
                            (1/2)*sigma^2* ((R[row,col+1,n_time]-2*R[row,col,n_time]+R[row,col-1,n_time])/(dv^2)) +
                            rho*sigma*sqrt(dv)*sqrt(col)*row*((R[row+1,col+1,n_time]-R[row-1,col+1,n_time]-R[row+1,col-1,n_time]+R[row-1,col-1,n_time])/(4*dv))) + 
                            R[row,col,n_time]
    }
  }
}


# Plot Results
# ---------------------------------------

levelpersp <- function(x, y, z, colors=heat.colors, ...) {
## getting the value of the midpoint
zz <- (z[-1,-1] + z[-1,-ncol(z)] + z[-nrow(z),-1] + z[-nrow(z),-ncol(z)])/4
## calculating the breaks
breaks <- hist(zz, plot=FALSE)$breaks
## cutting up zz
cols <- rev(colors(length(breaks)-1))
zzz <- cut(zz, breaks=breaks, labels=cols)
## plotting
persp(x, y, z, col=as.character(zzz), ...)
## return breaks and colors for the legend
list(breaks=breaks, colors=cols)
}

