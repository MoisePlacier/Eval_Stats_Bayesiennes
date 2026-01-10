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

data <- list( c= data$c, m= data$m , rm = data$rm, n = length(data$c) , Q= data$Q)


# ----------------------------------------------------------------------------
#                                    Model
# ----------------------------------------------------------------------------

# Read model file

source("model/model_hierarchical.txt")

# Write the model specification into a virtual text file, to be used by JAGS

model <- textConnection(model_string)


# ----------------------------------------------------------------------------
#                                  load inits
# ----------------------------------------------------------------------------

# Load inits values for 3 chains 

inits1 <- list(logit_theta = rep(0.1, data$n) , Nunif = rep(20000,data$n))

inits2 <- list( logit_theta = rep(0.5, data$n) , Nunif = rep(30000,data$n))

inits3 <- list( logit_theta = rep(0.9, data$n) , Nunif = rep(5000,data$n))

inits <- list(inits1,inits2,inits3)

# on initialise pas tau_z ni mu_z car JAGS s'en charge ! 


# ----------------------------------------------------------------------------
#                Set the vector for the variables to store
# ----------------------------------------------------------------------------

variables_set <- c("theta", "N", "logit_theta","a","b","sigma")


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
library(coda)

params <- c("a", "b", "sigma")

# Calcul du Gelman-Rubin diagnostic
gelman_res <- gelman.diag(mcmc[, params], confidence=0.95, transform=TRUE)

# Extraire les résultats
psrf <- gelman_res$psrf
df <- data.frame(
  Paramètre = rownames(psrf),
  PSRF_point = round(psrf[,1], 2),
  PSRF_upper95 = round(psrf[,2], 2),
  Interprétation = "Convergence satisfaisante"
)

# Ajouter ligne pour le multivarié
df <- rbind(df, data.frame(
  Paramètre = "Multivariate",
  PSRF_point = round(gelman_res$mpsrf, 2),
  PSRF_upper95 = NA,
  Interprétation = "Convergence globale"
))

# Afficher la table
print(df)
gelman.diag(mcmc[,c("a","b","sigma")], confidence=0.95, transform=TRUE)

# Pour quelques années représentatives de theta
gelman.diag(mcmc[,c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]","theta[6]","theta[7]","theta[8]","theta[9]","theta[10]","theta[11]","theta[12]","theta[13]","theta[14]","theta[15]","theta[16]","theta[17]","theta[18]")], confidence=0.95, transform=TRUE)

# Pour N
gelman.diag(mcmc[,c("N[1]","N[2]","N[3]","N[4]","N[5]","N[6]","N[7]","N[8]","N[9]","N[10]","N[11]","N[12]","N[13]","N[14]","N[15]","N[16]","N[17]","N[18]")], confidence=0.95, transform=TRUE)


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

# Extraire les échantillons MCMC pour theta et les hyperparamètres a, b, sigma
theta_logit <- as.matrix(mcmc[, grep("logit_theta", varnames(mcmc))])
a_samples <- as.vector(as.matrix(mcmc[,"a"]))
b_samples <- as.vector(as.matrix(mcmc[,"b"]))
sigma_samples <- as.vector(as.matrix(mcmc[,"sigma"]))

# Centrer et standardiser Q si nécessaire
Q <- data$Q
Q_scaled <- (Q - mean(Q)) / sd(Q)

# Nombre de chaînes et itérations
n_iter <- nrow(theta_logit)

# Matrice pour stocker les R2 pour chaque itération
R2_logit <- numeric(n_iter)

for(i in 1:n_iter){
  # Calcul de logit(theta) prédit par le modèle linéaire
  logit_theta_hat <- a_samples[i] + b_samples[i] * Q_scaled
  
  # Variance expliquée par le modèle (variance de logit_theta_hat)
  var_model <- var(logit_theta_hat)
  
  # Variance totale = variance expliquée + variance résiduelle sigma^2
  var_total <- var_model + sigma_samples[i]^2
  
  # R2-type
  R2_logit[i] <- var_model / var_total
}

# Résumé de R2
summary(R2_logit)
quantile(R2_logit, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))



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


library(ggplot2)

# Convertir l'objet mcmc.list en data.frame
mcmc.table <- as.data.frame(as.matrix(mcmc))

# Vérifier les noms des colonnes
head(colnames(mcmc.table))

# Extraire les colonnes correspondant à N[i]
N_cols <- grep("^N\\[", colnames(mcmc.table))

# Créer le tableau résumé
N_summary <- data.frame(
  year  = 1:length(N_cols),
  median = apply(mcmc.table[, N_cols], 2, median),
  lower  = apply(mcmc.table[, N_cols], 2, quantile, 0.05),
  upper  = apply(mcmc.table[, N_cols], 2, quantile, 0.95)
)

missing_years <- c(10,11)  # indices ou années exactes selon N_summary$year

library(ggplot2)

# Supposons que N_summary contient : year, median, lower, upper
# Et que missing_years est un vecteur des années manquantes, par exemple : 
missing_years <- c(10,11)  # indices ou années exactes selon N_summary$year

# Ajouter une colonne de couleur selon les années manquantes
N_summary$color <- ifelse(N_summary$year %in% missing_years, "red", "blue")

# Créer un plot avec médiane + intervalle crédible 5-95%
ggplot(N_summary, aes(x=factor(year))) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=color), width=0.2, size=1) +
  geom_point(aes(y=median, color=color), size=3) +
  scale_color_identity() +
  labs(
    x="Année", 
    y="N(t)", 
    title="Estimation du nombre de smolts migrants vers l'océan (IC : 95%)"
  ) +
  theme_minimal(base_size=14) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
