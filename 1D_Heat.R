# ================================================== 
# 1D Heat Equation: du/dt = k * (d**2u)/(dx**2)
# Model: M3
# Date: April 27, 2013
# Author: Jacob Roth (jroth13@cmc.edu)
# ==================================================



# INITIALIZING
# ------------------------------------------------- #
# Initialize function for temperature distribution  #
f <- function(l){                                   #
  return(sin(l))                                    #
}                                                   #
# ------------------------------------------------- #

# Set up the initial length and temperature of the thin rod 
l0 <- 0                                   # this is the x-coord of the low endpoint
ln <- 2*pi                                # this is the x-coord of the high endpoint
lsteps <- 50                              # specfy the number of intermediary steps between l0, ln
                                              # note this yields a vector of length lsteps + 1
dl <- (ln - l0) / lsteps                  # spacing between steps
l <- seq(l0, ln, by=dl)                   # generate discretized sequence representing l

temp0 <- f(l)                             # apply function f to specify initial temperature distribution (and time t)
temp0[1] <- templ0 <- 0                    # specify non-changing low endpoint through time
temp0[lsteps+1] <- templn <-1              # specify non-changing high endpoint through time

temp1 <- f(l)                             # generate vector to represent temperature distribution at time t +1
temp1 <- temp0


# Set up initial time steps 
t0 <- 0                                   # starting time for model analysis
tn <- 10000                                   # end time for model analysis
tsteps <- 10000                             # specify discretization of time into number of intermediary steps
dt <- (tn - t0) / tsteps                  # spacing between steps
t <- seq(t0, tn, by=dt)                   # discretized time sequence

# Make sure lambda <= 1
k <- .125                                   # parameter k in heat equation
lambda <- k * tsteps / lsteps^2           # see notes for derivation by discretization
print(lambda)          

# Build a matrix to track the resulting temperature across time
R <- matrix(NA, ncol=tsteps+1, nrow=lsteps+1)  # Results matrix will hold temperature for each odrered pair (l,t)


# RUN MODEL
# ===============================

R[,1] = temp0

# Calculate how temperature changes through time; temp0: at time t, temp1: at time t+1 
for (t in 2:tsteps){
  for (l in 2:lsteps){    
    temp1[l] <- temp0[l] + lambda * (temp0[l+1] + temp0[l-1] - 2 * temp0[l])  # along spatial point, we are... 
                                                                        # ... calcing temp at time t+1 from time t                         
  }
  
  for (l in 1:lsteps){
    temp0[l] = temp1[l]  # replace temp0 with temp1 to in preparation to run next iteration
  }
  R[,t] <- temp0     # once loop has run through each l at time n, store  in colun t + 1 of matrix R
}
print(R)

# ========================================
# PLOT RESULTS

plot(R[,1],type="l")
lines(R[,50])
lines(R[,100])
lines(R[,500])
lines(R[,1000])
lines(R[,500])

