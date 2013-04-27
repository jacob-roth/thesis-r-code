# ======================================= 
# Model: Explicit Black-Scholes 
# Date: April 27, 2013
# Author: Jacob Roth (jroth13@cmc.edu)
# Adapted from Mark Richardson's example
# =======================================

# Inputs
K = 75
r_yr = .05       # countinuously compounded rate
r = r_yr #r_yr/N # countinuously compounded rate
sigma_yr = .2    # yearly volatility
sigma = sigma_yr #sigma_yr * sqrt(dt)

M = 80           # asset nodes
N = 200          # time nodes
S_0 = 0          # minimum asset price
S_max = 150      # maximum asset price
t_min = 0        # minimum time
t_max = capT = 3 # maximum time (time at maturity)


# Variable Discretizations
S <- seq(S_0,S_max,by=S_max/M)   # asset discretization 
dS <- S_max/M

t <- seq(t_min,t_max,by=capT/N)  # backwards in time
tau = capT - t
dt <- capT/N

# ----------------------------
R = array(0,dim=c(M+1,N+1))

R[1,] = 0                              # BC at S_0
R[M+1,] = S_max - K * exp(-r*(tau))    # BC at S_max
R[,N+1] = pmax(S-K,0)                  # BC at maturity

# ----------------------------

# Define Coefficients A, B, and C:
A <- function(i){
  A_i <- (1/2) * dt * (r*i - sigma^2 * i^2)
  return(A_i)
}

B <- function(i){
  B_i <- 1 + dt * (sigma^2 * i^2 + r)
  return(B_i)
}

C <- function(i){
  C_i <- (1/2) * dt * (-r*i - sigma^2 * i^2)
  return(C_i)
}

# ----------------------------

Coeff = array(0,dim=c(M+1,M+1))  # initialize coeff matrix, dim (M+1)x(M+1) to match (M+1) stock soln vector
a = rep(0,times=M+1)
b = rep(0,times=M+1)                   
c = rep(0,times=M+1)
# a = A(seq(1,M))   # sub-diagonal, a[col]
# b = B(seq(0,M))   # main-diagonal
# c = C(seq(-1,M))  # super-diagonal, c[col]

for (i in 1:(M+1)){ # build vector of coeffs; evaluate at i-1 to account for next loop / matrix starting at 1 
  a[i] = A(i-1)
  b[i] = B(i-1)
  c[i] = C(i-1)  
}

diag(Coeff) = b     # input b onto main diagonal

for (col in 1:(M+1)){ # input 'a' coeff vector onto lower band
  for(row in 1:(M+1)){
    if (row == col+1){
      Coeff[row,col] = a[row]
    }
  }
}

for (col in 1:(M+1)){ # input 'c' coeff vector onto upper band
  for(row in 1:(M+1)){
    if (row == col-1){
      Coeff[row,col] = c[row]
    }
  }
}

Coeff_inv = solve(Coeff) # invert coeff matrix

for (col in N:1){               # for each new time step; start with t_max - 1 (i.e. N) to work off BC at N+1 
  S_u = rep(0,times=M+1)        # create vector to hold S_u, the unknown values of S at each time step from N:1
  S_u[1] = A(0) * R[1,col+1] # prepare S_u[1], i.e. the A_0 that is known and gets moved over to the RHS
  S_u[length(S_u)] = C(M+1) * R[M+1,col+1]  # similarly, prepare to move known quantity S_u[end] over to RHS 
  RHS =  R[,col+1] - S_u     # we're solving for RHS which we found in last step, so RHS_new = RHS_old with S_u adjustment
  S_u = Coeff_inv %*% RHS       # solve
  R[(2:M),col] = S_u[2:M]    # and store in appropriate column
}

shading <- rev(palette(rainbow(100)))
jet.colors <- colorRampPalette(shading[28:100])

nbcol <- 1000
color <- jet.colors(nbcol)
zfacet <- R[-1,-1] #+ R[-1,-(M+1)] + R[-(N+1),-1] + R[-(N+1),-1]
facetcol <- cut(zfacet,nbcol)
persp(R,col=color[facetcol],phi=10,theta=-50,shade=0.001,ltheta=100,
      xlab="Stock Price",ylab="Tau",zlab="Option Value",ticktype="detailed")
require(rgl)
persp3d(R,col=color['red'],phi=10,theta=-28,shade=0.001,ltheta=100,
        xlab="Stock Price",ylab="Tau",zlab="Option Value",ticktype="detailed")
