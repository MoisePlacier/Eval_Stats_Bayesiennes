
rm(list = ls())

# Load library

library(rjags)
library(coda)

# load module DIC

load.module("dic")

data <- read.table("data/data_CMR_smolts.txt", header = TRUE)
data <- data[!is.na(data$Nb_juv0_previous_year), ]
data

# Build data set to be read in JAGS


data <- list( c= data$c, m= data$m , rm = data$rm, n = length(data$c) , Q= data$Q, pSm1 = data$pSm1, Njuv_obs = data$Nb_juv0_previous_year)

# Read model file

source("model/model_survival_rate.txt")

model <- textConnection(model_string)

# Load inits values for 3 chains 

inits1 <- list(logit_theta = rep(0.1, data$n), Nunif = rep(20000, data$n), surv = rep(0.5, data$n))
inits2 <- list(logit_theta = rep(0.5, data$n), Nunif = rep(30000, data$n), surv = rep(0.6, data$n))
inits3 <- list(logit_theta = rep(0.9, data$n), Nunif = rep(5000, data$n), surv = rep(0.4, data$n))

inits <- list(inits1,inits2,inits3)


# ----------------------------------------------------------------------------
#                Set the vector for the variables to store
# ----------------------------------------------------------------------------

variables_set <- c("theta", "N", "logit_theta","a","b","sigma", "N_smolt1","surv","Njuv")


# ----------------------------------------------------------------------------
#                      COMPILE and ADAPT MCMC samplers
# ----------------------------------------------------------------------------

n_chains = 3
n_adapt = 10000 # "période de chauffe" avant de  

time_compile <- system.time(
model_jags <- jags.model(model, data=data, inits=inits, n.chain=n_chains, n.adapt=n_adapt)
)
time_compile


# ----------------------------------------------------------------------------
#                        RUN and STORE MCMC samples
# ----------------------------------------------------------------------------

thin = 10
n_samples = 1000
n_iter = n_samples*thin

# stockage de 1/10 (thin = 10) des valeurs de theta et N pour pas surcharger 

time_compile <- system.time(
mcmc <- coda.samples(model=model_jags, variable.names = variables_set, n.iter=n_iter, thin=thin)
)
time_compile



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

params <- c("surv[7]", "surv[10]", "surv[2]")

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
gelman.diag(mcmc[,c("surv[1]","surv[2]","surv[3]","surv[4]","surv[5]","surv[6]","surv[7]","surv[8]","surv[9]","surv[10]","surv[11]","surv[12]","surv[13]","surv[14]","surv[15]","surv[16]","surv[17]")], confidence=0.95, transform=TRUE)

# Pour N
gelman.diag(mcmc[,c("Njuv[1]","Njuv[2]","Njuv[3]","Njuv[4]","Njuv[5]","Njuv[6]","Njuv[7]","Njuv[8]","Njuv[9]","Njuv[10]","Njuv[11]","Njuv[12]","Njuv[13]","Njuv[14]","Njuv[15]","Njuv[16]","Njuv[17]")], confidence=0.95, transform=TRUE)


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
missing_years <- c(9,10)  # indices ou années exactes selon N_summary$year

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





# Convertir l'objet mcmc.list en data.frame
mcmc.table <- as.data.frame(as.matrix(mcmc))

# Vérifier les noms des colonnes
colnames(mcmc.table)

# Extraire les colonnes correspondant à N[i]
N_cols <- grep("^surv\\[", colnames(mcmc.table))

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
missing_years <- c(9,10)  # indices ou années exactes selon N_summary$year

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
library(ggplot2)
library(dplyr)
library(tidyr)

# Convertir l'objet mcmc.list en data.frame
mcmc.table <- as.data.frame(as.matrix(mcmc))

# Conversion en format long
N_long <- data.frame(
  year = rep(1:length(N_cols), each = nrow(mcmc.table)),
  value = as.vector(as.matrix(mcmc.table[, N_cols])),
  type = "N_total"
)

Njuv_long <- data.frame(
  year = rep(1:length(Njuv_cols), each = nrow(mcmc.table)),
  value = as.vector(as.matrix(mcmc.table[, Njuv_cols])),
  type = "N_juv"
)

df_long <- bind_rows(N_long, Njuv_long)

# Boxplot
ggplot(df_long, aes(x = factor(year), y = value, fill = type, color = type)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +  # remplissage semi-transparent, cacher outliers pour plus de clarté
  scale_fill_manual(values = c("N_total" = "blue", "N_juv" = "orange")) +
  scale_color_manual(values = c("N_total" = "darkblue", "N_juv" = "darkorange")) + # contour plus foncé
  labs(
    x = "Année",
    y = "Abondance estimée",
    title = "Distribution annuelle du nombre de juvéniles et de smolts",
    fill = "Type",
    color = "Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



library(ggplot2)
library(dplyr)
library(tidyr)

# Supposons df_long contient N_total, N_juv et surv

# Conversion MCMC en data.frame
mcmc.table <- as.data.frame(as.matrix(mcmc))

# Colonnes Njuv et surv
Njuv_cols <- grep("^Njuv\\[", colnames(mcmc.table))
Nsmolt1_cols <- grep("^N_smolt1\\[", colnames(mcmc.table))
surv_cols <- grep("^surv\\[", colnames(mcmc.table))
N_cols <- grep("^N\\[", colnames(mcmc.table))

# Passage en format long
dd_long <- data.frame(
  year = factor(rep(1:length(Njuv_cols), each = nrow(mcmc.table))),
  Njuv = as.vector(as.matrix(mcmc.table[, Njuv_cols])),
  N = as.vector(as.matrix(mcmc.table[, N_cols])),
  surv = as.vector(as.matrix(mcmc.table[, surv_cols])),
  Nsmolt1 = as.vector(as.matrix(mcmc.table[, Nsmolt1_cols]))
)

N_long <- data.frame(
  year = rep(1:length(N_cols), each = nrow(mcmc.table)),
  value = as.vector(as.matrix(mcmc.table[, N_cols])),
  type = "N"
)

Njuv_long <- data.frame(
  year = rep(1:length(Njuv_cols), each = nrow(mcmc.table)),
  value = as.vector(as.matrix(mcmc.table[, Njuv_cols])),
  type = "N_juv"
)

Nsmolt1_long <- data.frame(
  year = rep(1:length(Njuv_cols), each = nrow(mcmc.table)),
  value = as.vector(as.matrix(mcmc.table[, Nsmolt1_cols])),
  type = "Nsmolt1"
)

df_long <- bind_rows(N_long, Njuv_long,Nsmolt1_long)


# Transformer en format long pour ggplot
df_plot <- df_plot %>%
  select(year, median, lower, upper, type) %>%
  mutate(year = factor(year))


# Pour superposer sur la même figure, il faut mettre la survie sur une échelle comparable
# On normalise ici entre 0 et max(median N_total/N_juv)
max_abund <- max(df_plot$upper)

# Préparer les données de survie
surv_summary <- dd_long %>%
  group_by(year) %>%
  summarise(
    median = median(surv),
    lower = quantile(surv, 0.05),
    upper = quantile(surv, 0.95)
  ) %>%
  mutate(
    year_num = as.numeric(as.character(year)),   # X numérique pour geom_ribbon
    median_scaled = median * max_abund,
    lower_scaled = lower * max_abund,
    upper_scaled = upper * max_abund
  )

# Plot combiné
ggplot() +
  # Boxplots N_total et N_juv
  geom_boxplot(data = df_long %>% filter(type %in% c("N","N_juv","Nsmolt1")),
               aes(x = factor(year), y = value, fill = type, color = type),
               alpha = 0.4, outlier.shape = NA) +
  scale_fill_manual(values = c("N" = "blue","Nsmolt1"="green", "N_juv" = "orange")) +
  scale_color_manual(values = c("N" = "darkblue", "N_juv" = "darkorange","Nsmolt1"="darkgreen")) +
  
  # Survie en ligne + intervalle (ribbon)
  geom_ribbon(data = surv_summary,
              aes(x = year_num, ymin = lower_scaled, ymax = upper_scaled),
              fill = "red", alpha = 0.2) +
  geom_line(data = surv_summary,
            aes(x = year_num, y = median_scaled),
            color = "red", size = 1) +
  
  # Axes
  scale_y_continuous(
    name = "Abondance estimée",
    sec.axis = sec_axis(~ ./max_abund, name = "Taux de survie")
  ) +
  
  labs(
    x = "Année",
    title = "Abondance annuelle et distribution de la survie"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######################################################
############ distrib surv
########################################
library(ggplot2)

surv_long <- data.frame(
  year = factor(rep(seq_along(surv_cols), each = nrow(mcmc.table))),
  surv = as.vector(as.matrix(mcmc.table[, surv_cols]))
)

# Années avec données juvéniles manquantes
missing_years <- c(9, 10)

# Indicateur logique (plus robuste que des couleurs codées à la main)
surv_long$missing <- surv_long$year %in% missing_years

# -----------------------
# Boxplots
# -----------------------
ggplot(surv_long, aes(x = year, y = surv, fill = missing)) +
  geom_boxplot(outlier.size = 0.6) +
  scale_fill_manual(
    values = c("FALSE" = "grey80", "TRUE" = "red"),
    labels = c("FALSE" = "Données complètes", "TRUE" = "données rm manquantes")
  ) +
  labs(
    x = "Année",
    y = "Taux de survie juvéniles → smolts",
    fill = "",
    title = "Distribution a posteriori du taux de survie annuel"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
# -----------------------
# Variante 2 : Violin plots
# -----------------------
ggplot(surv_long, aes(x = year, y = surv)) +
  geom_violin(fill = "grey80", color = "black", trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  labs(
    x = "Année",
    y = "Taux de survie juvéniles → smolts",
    title = "Distribution a posteriori du taux de survie annuel"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



library(ggplot2)

# Conversion MCMC en data.frame
mcmc.table <- as.data.frame(as.matrix(mcmc))

# Colonnes Njuv et surv
Njuv_cols <- grep("^Njuv\\[", colnames(mcmc.table))
surv_cols <- grep("^surv\\[", colnames(mcmc.table))

# Passage en format long
dd_long <- data.frame(
  year = factor(rep(1:length(Njuv_cols), each = nrow(mcmc.table))),
  log_Njuv = log(as.vector(as.matrix(mcmc.table[, Njuv_cols]))),
  surv = as.vector(as.matrix(mcmc.table[, surv_cols]))
)

ggplot(dd_long, aes(x = log_Njuv, y = surv)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ year, scales = "free_x") +
  labs(
    x = expression(log(N[juv](t))),
    y = "surv(t)",
    title = "Relation densité–survie par année (a posteriori)"
  ) +
  theme_minimal(base_size = 12)
library(dplyr)

cor_post <- dd_long %>%
  group_by(year) %>%
  summarise(
    cor = cor(log_Njuv, surv)
  )

cor_post <- dd_long %>%
  group_by(year) %>%
  summarise(
    cor_median = median(sapply(1:n(), function(i) cor(log_Njuv[i], surv[i]))),  # exemple, à remplacer par vrai bootstrap
    cor_lower  = quantile(sapply(1:n(), function(i) cor(log_Njuv[i], surv[i])), 0.05),
    cor_upper  = quantile(sapply(1:n(), function(i) cor(log_Njuv[i], surv[i])), 0.95)
  )

ggplot(cor_post, aes(x = year, y = cor)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Année",
    y = "Corrélation postérieure",
    title = "Corrélation a posteriori entre log(Njuv) et surv"
  ) +
  theme_minimal(base_size = 14)



library(ggplot2)
library(dplyr)

# Colonnes MCMC pour Njuv et N_smolt1
Njuv_cols <- grep("^Njuv\\[", colnames(mcmc.table), value = TRUE)
Nsmolt1_cols <- grep("^N_smolt1\\[", colnames(mcmc.table), value = TRUE)
surv_cols <- grep("^surv\\[", colnames(mcmc.table), value = TRUE)

# Exemple pour les années 1 à 5 (ou toutes)
years <- c(1,2,3,4,7,8,9,10,14)

# Construire un data.frame "long" avec Njuv, Nsmolt1 et surv pour chaque année
scatter_df <- bind_rows(lapply(years, function(i) {
  data.frame(
    year = i,
    Njuv = mcmc.table[[paste0("Njuv[",i,"]")]],
    N_smolt1 = mcmc.table[[paste0("N_smolt1[",i,"]")]],
    surv = mcmc.table[[paste0("surv[",i,"]")]]
  )
}))

# Calcul de la survie médiane par année
surv_median <- scatter_df %>%
  group_by(year) %>%
  summarise(surv_med = median(surv))

# Plot
ggplot(scatter_df, aes(x = N_smolt1, y = Njuv)) +
  geom_point(aes(color = surv), alpha = 0.4) +
  geom_abline(data = surv_median, aes(slope = 1/surv_med, intercept = 0),
              color = "red", linetype = "dashed", size = 1) +
  facet_wrap(~year, scales = "free") +
  scale_color_viridis_c(option = "plasma") +
  labs(
    x = "Nombre de smolts d'un an (N_smolt1)",
    y = "Nombre de juvéniles latents estimés(Njuv)",
    color = "Survie medianne (postérieure)",
    title = "Relation postérieure Njuv ↔ N_smolt1 par année",
    subtitle = "Chaque point = un échantillon MCMC ; ligne rouge = survie médiane"
  ) +
  theme_minimal(base_size = 14)

ggplot(dd_long, aes(x = log_Njuv, y = surv)) +
  geom_point(alpha = 0.3, color = "blue") +  # points transparents pour visualiser la densité
  geom_smooth(method = "loess", color = "red") +  # ligne de tendance locale
  labs(x = "log(Njuv)", y = "Taux de survie", 
       title = "Relation entre log(Njuv) et taux de survie") +
  theme_minimal(base_size = 14)
ggplot(dd_long, aes(x = log_Njuv, y = surv, color = year)) +
  geom_point(alpha = 0.5) +             # points semi-transparents pour visualiser la densité
  geom_smooth(method = "loess", color = "black", se = FALSE) +  # ligne de tendance locale globale
  scale_color_viridis_d(name = "Année") +  # palette discrète pour les années
  labs(x = "log(Njuv)", y = "Taux de survie", 
       title = "Relation entre log(Njuv) et taux de survie par année") +
  theme_minimal(base_size = 14)

dd_long$flow <- rep(as.factor(data$Q), each = nrow(mcmc.table))


# Boxplot
ggplot(dd_long, aes(x = flow, y = surv)) +
  geom_boxplot(alpha = 0.5, fill = "skyblue", color = "blue") +
  labs(x = "Débit", y = "Taux de survie",
       title = "Relation entre débit et taux de survie") +
  theme_minimal(base_size = 14)
