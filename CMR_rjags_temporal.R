# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Model for smolt run estimation
# E. Rivot - November 2025
# ---------------------------------------------------------------
# ---------------------------------------------------------------


# ---------------------------------------------------------------

# Check the path to JAGS.exe

# Optional if rjags can not connect directly to the JAGS.exe installed on your computer
# You will see that in the information message after loading the rjags library
# CAUTION : the path C:\\ is to be adapted to redirect where your JAGS.exe is installed
# the name JAGS-4.2.0 should also be changed if another version of JAGS is installed on your computer

# Sys.setenv(JAGS_HOME = "C:\\prog_td\\JAGS\\JAGS-4.2.0\\")

# ---------------------------------------------------------------

rm(list = ls())

# Load library

library(rjags)
library(coda)

# load module DIC

load.module("dic")


# ----------------------------------------------------------------------------
#                                 Load data
# ----------------------------------------------------------------------------

data <- read.table("data/data_CMR_smolts.txt", header = TRUE)
data

# Build data set to be read in JAGS
# -------------------------------

data <- list( c= data$c, m= data$m , rm = data$rm, n = length(data$c) )


# ----------------------------------------------------------------------------
#                                    Model
# ----------------------------------------------------------------------------

# Read model file

source("model/model_temporal.txt")

# Write the model specification into a virtual text file, to be used by JAGS

model <- textConnection(model_string)


# ----------------------------------------------------------------------------
#                                  load inits
# ----------------------------------------------------------------------------

# Load inits values for 3 chains 

inits1 <- list(theta = rep(0.1, data$n) , Nunif = rep(20000,data$n))

inits2 <- list( theta = rep(0.5, data$n) , Nunif = rep(30000,data$n))

inits3 <- list( theta = rep(0.9, data$n) , Nunif = rep(5000,data$n))

inits <- list(inits1,inits2,inits3)


# ----------------------------------------------------------------------------
#                Set the vector for the variables to store
# ----------------------------------------------------------------------------

variables_set <- c("theta", "N")


# ----------------------------------------------------------------------------
#                      COMPILE and ADAPT MCMC samplers
# ----------------------------------------------------------------------------
# Compile the model for n_chains
# and adapt the MCMC samplers during n.adapt iterations
# + record the computer time needed to run the model

n_chains = 3
n_adapt = 10000 # "période de chauffe" avant de  

# Warning adapting - if rjags return a warning "adapting incomplete", this means
# that MCMC sampler will not be optimized
# --> need to increase n.adapt so as the warning message desapears

time_compile <- system.time(
model_jags <- jags.model(model, data=data, inits=inits, n.chain=n_chains, n.adapt=n_adapt)
)
time_compile


# ----------------------------------------------------------------------------
#                        RUN and STORE MCMC samples
# ----------------------------------------------------------------------------
# Run the model for a number of iteration = samples
# thin the chain (keep only one out of "thin" iterations
# record the computer time needed to run the model

thin = 10
n_samples = 1000
n_iter = n_samples*thin

# stockage de 1/10 (thin = 10) des valeurs de theta et N pour pas surcharger 

time_compile <- system.time(
mcmc <- coda.samples(model=model_jags, variable.names = variables_set, n.iter=n_iter, thin=thin)
)
time_compile


# Calculate DIC with n.DIC new (additional) MCMC iterations
# n_DIC = 1000
# dic.samples(model = model_jags, n.iter = n_DIC, thin = 1, type = "pD")




# ----------------------------------------------------------------------------
#                             Explore the results
# ----------------------------------------------------------------------------

# --------------------------------------------------------
# Work with mcmc.list
# --------------------------------------------------------

# "mcmc" is an object of the class "mcmc.list" (see package library(coda)
# to explore, plot ... mcmc objects

is(mcmc)

# Names and dimention of the variables stored in the mcmc list

varnames(mcmc)
nvar(mcmc)

nchain(mcmc)
niter(mcmc)


# ---------------------------------------------------
# Example of code to extract the variable of interest from the mcmc list
# ---------------------------------------------------

# mcmc.var is still an mcmc list but contains only the variable of interest

var = "theta"
mcmc.var = mcmc[,which(substr(varnames(mcmc),1,nchar(var))==paste(var,"",sep=""))]

# Extract iterations from start = 100 to end = 200
# window() is a function of the coda package design to manage mcmc.list object


# --------------------------------------------------------
# Assess convergence and effective sample size
# --------------------------------------------------------

# Gelman-Rubin convergence diagnostics
# Point est. should be near 1


# Gelman-Rubin convergence diagnostics
# Point est. should be near 1
# TryCatch will avoid stoping the program whn an error occurs in gelman.diag function 
tryCatch(eval(
  
  gelman.diag(mcmc, confidence = 0.95, transform=TRUE, autoburnin=TRUE)
  # You may alternatively Select some particular nodes to test
  # gelman.diag(mcmc[,c("theta","N")], confidence = 0.95, transform=TRUE, autoburnin=TRUE)
  # gelman.diag(mcmc[,varnames(mcmc)[1:nvar(mcmc)]], confidence = 0.95, transform=TRUE, autoburnin=TRUE)
  
), 
error = function(e) message("Oops!  ", as.character(e)))


# Effective sample size - should ideally be >1000 for all variables
# TryCatch will avoid stoping the program when an error occurs in effectiveSize function
tryCatch(eval(
  
  {print("Relative_ESS") ; effectiveSize(mcmc)/(n_samples*n_chains)}
  
), 
error = function(e) message("Oops!  ", as.character(e)))



# --------------------------------------------
# Get summary of estimates
# --------------------------------------------

# Summary statistics
tryCatch(eval(
  
  summary(mcmc, quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95))[2]
  
), 
error = function(e) message("Oops!  ", as.character(e)))



# --------------------------------------------
# plot results
# --------------------------------------------

# source("plot_results.R")

# Use this code to run the plot_results program with no errors
plot_results <- parse(file = "plot_results.R")
for (i in seq_along(plot_results)) {
  tryCatch(eval(plot_results[[i]]), 
           error = function(e) message("Oops!  ", as.character(e)))
}






