# ===================================== 
# German's Stochastic Volatility Model 
# Model: M3
# Date: April 27, 2013
# Author: Jacob Roth (jroth13@cmc.edu)
# =====================================


# Discretize Variables
# -------------------------------------
#   Asset:
S_0 <-   0                                   
S_max <- 150                                   
Ssteps <- L <- 10                              
S <- l <- seq(S_0, S_max, length.out=Ssteps)   
dS = S[2] - S[1]

#   Volatility:
v_0 <- 0
v_max <- 1
vsteps <- W <- 12
v <- w <- seq(v_0, v_max, length.out=vsteps)
dv <- v[2] - v[1]

#   Time:
t_0 <- 0
t_max <- capT <- .01
tsteps <- N <- 100
t <- n <- seq(t_max, t_0, length.out=tsteps)
tau = capT- t
dt <- t_max / (tsteps-1)

#   Model Parameters
K1 <- 75
K2 <- 125
rho <- .1
r <- .1
b <- .1
a <- 0.1
mu <- .3
sigma <- 3
Vsum <- 0

# Discretization Functions
V1L <- function(R,x,y,t,mu){
  return(mu*(R[x,y+1,t] - R[x,y,t])/dv)
}
V1R <- function(R,x,y,t,mu){
  return(mu*(R[x,y,t] - R[x,y-1,t])/dv)
}
V1M <- function(R,x,y,t,mu){
  return(mu*(R[x,y+1,t] - R[x,y-1,t])/(2*dv))
}
# ------------------------
S2L <- function(R,x,y,t){
  return(.5*y*dv*(x*dS)^2*(R[x+2,y,t] - 2*R[x+1,y,t] + R[x,y,t])/(dS)^2)
}
S2R <- function(R,x,y,t){
  return(.5*y*dv*(x*dS)^2*(R[x,y,t] - 2*R[x-1,y,t] + R[x-2,y,t])/(dS)^2)
}
S2M <- function(R,x,y,t){
  return(.5*y*dv*(x*dS)^2*(R[x+1,y,t] - 2*R[x,y,t] + R[x-1,y,t])/(dS)^2)
}
# ------------------------
V2L <- function(R,x,y,t,sigma){
  return(.5*sigma*(R[x,y+2,t] - R[x,y+1,t] + R[x,y,t])/dv^2)
}
V2R <- function(R,x,y,t,sigma){
  return(.5*sigma*(R[x,y,t] - R[x,y-1,t] + R[x,y-2,t])/dv^2)
}
V2M <- function(R,x,y,t,sigma){
  return(.5*sigma*(R[x,y+1,t] - R[x,y,t] + R[x,y-1,t])/dv^2)
} 
# ------------------------
SVL <- function(R,x,y,t,sigma,rho){
  return(rho*(sigma*sqrt(y*dv)*(x*dS))*(R[x+1,y+1,t] - R[x+1,y,t] - R[x,y+1,t] + R[x,y,t])/(dS*dv))
}
SVLxRy <- function(R,x,y,t,sigma,rho){
  return(rho*(sigma*sqrt(y*dv)*(x*dS))*(R[x+1,y,t] - R[x+1,y-1,t] - R[x,y,t] + R[x,y-1,t])/(dS*dv))
}
SVRxLy <- function(R,x,y,t,sigma,rho){
  return(rho*(sigma*sqrt(y*dv)*(x*dS))*(R[x,y+1,t] - R[x,y,t] - R[x-1,y+1,t] + R[x-1,y,t])/(dS*dv))
}
SVR <- function(R,x,y,t,sigma,rho){
  return((rho*sigma*sqrt(y*dv)*(x*dS))*(R[x-1,y-1,t] - R[x-1,y,t] - R[x,y-1,t] + R[x,y,t])/(dS*dv))
}
SVM <- function(R,x,y,t,sigma,rho){
  return(rho*(sigma*sqrt(y*dv)*(x*dS))*(R[x+1,y+1,t] - R[x+1,y-1,t] - R[x-1,y+1,t] + R[x-1,y-1,t])/(4*dS*dv))
}

# Initialize Boundary Conditions in Solution Arrays
#   Stock/Volatility Grid:
Uold <- array(NA,dim=c(L,W),dimnames=c("asset","vol"))
Unew <- array(NA,dim=c(L,W),dimnames=c("asset","vol"))
R <- array(0,dim=c(L,W,N))

#   Payoff at Maturity (t = capT):
for (x in 1:L){
  for (y in 1:W){
    Uold[x,y] <- pmin(pmax(S[x] - K1, 0),K2 - K1) 
  }
}
#   t = capT BC
R[,,N] = Uold                                                     

# Mesh Calculation
# --------------------------------
for(n_time in (N:2)){
  for(col in (W:1)){
    for(row in (L:1)){
      location <- c(row,col)
      
      Vsum <- Vsum + col*dv
      b <- Vsum/n_time
      mu <- a*(b-dv*col)
      sigma <- sigma * sqrt(dv*col)
      
      if(location[1] <= 2){
        if(location[2] >= (W-2)){
          formula <- 1
        }
        if(location[2] > 2 && location[2] < (W-2)){
          formula <- 2
        }
        if(location[2] <= 2){
          formula <- 3
        }
      }
      if(location[1] > 2 && location[1] < (W-2)){
        if(location[2] >= (W-2)){
          formula <- 4
        }
        if(location[2] < (W-2) && location[2] > 2){
          formula <- 5
        }
        if(location[2] <= 2){
          formula <- 6
        }
      }
      if(location[1] >= (L-2)){
        if(location[2] >= (W-2)){
          formula <- 7
        }
        if(location[2] < (W-2) && location[2] > 2){
          formula <- 8
        }
        if(location[2] <= 2){
          formula <- 9
        }
      }
      if(formula == 1){
        R[row,col,n_time-1] <- (dt*(V1R(R,row,col,n_time,mu) + S2L(R,row,col,n_time) + V2R(R,row,col,n_time,sigma) + SVLxRy(R,row,col,n_time,sigma,rho)) + R[row,col,n_time]) 
      }
      if(formula == 2){
        R[row,col,n_time-1] <- (dt*(V1M(R,row,col,n_time,mu) + S2L(R,row,col,n_time) + V2M(R,row,col,n_time,sigma) + SVLxRy(R,row,col,n_time,sigma,rho)) + R[row,col,n_time])
      }
      if(formula == 3){
        R[row,col,n_time-1] <- (dt*(V1L(R,row,col,n_time,mu) + S2L(R,row,col,n_time) + V2L(R,row,col,n_time,sigma) + SVL(R,row,col,n_time,sigma,rho)) + R[row,col,n_time])
      }
      if(formula == 4){
        R[row,col,n_time-1] <- (dt*(V1R(R,row,col,n_time,mu) + S2M(R,row,col,n_time) + V2R(R,row,col,n_time,sigma) + SVR(R,row,col,n_time,sigma,rho)) + R[row,col,n_time])
      }
      if(formula == 5){
        R[row,col,n_time-1] <- (dt*(V1M(R,row,col,n_time,mu) + S2M(R,row,col,n_time) + V2M(R,row,col,n_time,sigma) + SVM(R,row,col,n_time,sigma,rho)) + R[row,col,n_time])
      }
      if(formula == 6){
        R[row,col,n_time-1] <- (dt*(V1L(R,row,col,n_time,mu) + S2M(R,row,col,n_time) + V2L(R,row,col,n_time,sigma) + SVL(R,row,col,n_time,sigma,rho)) + R[row,col,n_time])
      }
      if(formula == 7){
        R[row,col,n_time-1] <- (dt*(V1R(R,row,col,n_time,mu) + S2R(R,row,col,n_time) + V2R(R,row,col,n_time,sigma) + SVR(R,row,col,n_time,sigma,rho)) + R[row,col,n_time])
      }
      if(formula == 8){
        R[row,col,n_time-1] <- (dt*(V1M(R,row,col,n_time,mu) + S2R(R,row,col,n_time) + V2M(R,row,col,n_time,sigma) + SVRxLy(R,row,col,n_time,sigma,rho)) + R[row,col,n_time])
      }
      if(formula == 9){
        R[row,col,n_time-1] <- (dt*(V1L(R,row,col,n_time,mu) + S2R(R,row,col,n_time) + V2L(R,row,col,n_time,sigma) + SVRxLy(R,row,col,n_time,sigma,rho)) + R[row,col,n_time])
      }
    }
  }
}

# Plotting Results
# ---------------------------------
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

levelpersp(S,v,R[,,100],xlab="Asset Price", theta=-30,
           ylab="Volatility",zlab="Option Value",ticktype="detailed", 
           main="Option Value Over Volatility",phi=30) # plots val across vols at time)
levelpersp(S[4:8],v[2:9],R[4:8,2:9,1],theta=-30,xlab="Asset Price",
           ylab="Time",zlab="Option Value",ticktype="detailed", 
           main="Zoom Option Value at t = 0.1")

R[5,4,1] # option value at 75, 30% vol -- change S_0 to 10
R[6,4,1] # option value at 75, 30% vol -- change S_0 to 0
 
