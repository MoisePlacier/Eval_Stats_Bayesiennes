# ---------------------------------------------------------------
# ---------------------------------------------------------------
# PLOT RESULTS
# Model for smolt run estimation
# E. Rivot - November 2025
# ---------------------------------------------------------------
# ---------------------------------------------------------------


# Source graphic functions
# ------------------------------------------------------------------------

source("function/f.panel.cor.R")
source("function/f.panel.dens.R")
source("function/f.density.bivar.R")


# ----------------------------------------------------------------------------
#                             Explore the results
# ----------------------------------------------------------------------------


# ---------------------------------------------------
# Extract the variable of interest from the mcmc list
# ---------------------------------------------------

# mcmc.var is still an mcmc list but contains only the variable of interest

var = "theta"
mcmc.var = mcmc[,which(substr(varnames(mcmc),1,nchar(var))==paste(var,"",sep=""))]

# Extract iterations from start = 100 to end = 200
# quartz() is a function of the coda package design to manage mcmc.list object

# start = 100
# end = 200
# mcmc <- quartz(mcmc, start=start,end=end)



# --------------------------------------------
# Work with mcmc samples stored in TABLES
# --------------------------------------------

# Extract MCMC chains and store in a TABLE
# sometimes easier to manipulate than mcmc.list object
# each column = all mcmc samples for each variable

# as.matrix() does not work if a multidimensional variable is in results
# here, "res" is multidimensionnal
# works after quartz is used

mcmc.table <- as.data.frame(as.matrix(quartz(mcmc)))
mcmc.table <- as.data.frame(as.matrix(mcmc))
head(mcmc.table)



# Traceplot of MCMC chains (see also ?traceplot - coda package)
# -------------------------------------------------------------

# Trace plot all variables (works wathever the dimension)
quartz()
plot(mcmc, trace = TRUE, density = FALSE)

# Other synthax - Will not work if the variables are vectors
# quartz()
# par(mfrow=c(2,1))
# traceplot(mcmc[,'theta'],ylab="theta")
# traceplot(mcmc[,'N'],ylab="N")


# Density plot
# -----------------------------------

# Density plot all variables (works wathever the dimension)
quartz()
plot(mcmc, trace = FALSE, density = TRUE)

# Other synthax - Will not work if the variables are vectors
# quartz()
# par(mfrow=c(2,1))
# densplot(mcmc[,'theta'],xlab="theta")
# densplot(mcmc[,'N'],xlab="N")


# Trace plot and density plot 
# -----------------------------------

# Theta[1] and N[1]
# Will only work if dim theta and N > 1

if (length(as.vector(mcmc.table$'theta[1]'))!=0){
  
quartz()
par(mfrow=c(2,2))
traceplot(mcmc[,'theta[1]'],ylab="theta")
traceplot(mcmc[,'N[1]'],ylab="N")

densplot(mcmc[,'theta[1]'],xlab="theta")
densplot(mcmc[,'N[1]'],xlab="N")
} else {
  quartz()
  par(mfrow=c(2,2))
  traceplot(mcmc[,'theta'],ylab="theta")
  traceplot(mcmc[,'N'],ylab="N")
  
  densplot(mcmc[,'theta'],xlab="theta")
  densplot(mcmc[,'N'],xlab="N")  
}


# Parameters a,b,sigma
# Will only work if a exists in the model

if (length(as.vector(mcmc.table$'a'))!=0){
  
quartz()
par(mfrow=c(2,3))
traceplot(mcmc[,'a'],ylab="a")
traceplot(mcmc[,'b'],ylab="b")
traceplot(mcmc[,'sigma'],ylab="sigma")

densplot(mcmc[,'a'],ylab="a")
densplot(mcmc[,'b'],ylab="b")
densplot(mcmc[,'sigma'],ylab="sigma")
}



# --------------------------------------------
# Loi jointe 2D de (theta,N) en une dimension
# --------------------------------------------

library(MASS) # necessaire pour la fonction kde2d qui permet de faire les contours

# Will only work if dim theta and N = 1
if (length(as.vector(mcmc.table$'theta'))!=0){
  dim1 <- as.vector(mcmc.table$'theta')
  dim2 <- as.vector(mcmc.table$'N')
} else {
  # Will only work if dim theta and N > 1
  dim1 <- as.vector(mcmc.table$'theta[1]')
  dim2 <- as.vector(mcmc.table$'N[1]')  
}

def.par <- par(no.readonly = TRUE)
quartz()
f.density.bivar(dim1,dim2,nlevels=3,nb.points=2000)


# --------------------------------------------
#              Boxplot theta and N
# --------------------------------------------

# Boxplot theta
# --------------------------------------------

# Will only work if dim theta and N > 1

if (length(as.vector(mcmc.table$'theta[1]'))!=0){
  
var = "theta"
mcmc.var = mcmc[,which(substr(varnames(mcmc),1,nchar(var))==paste(var,"",sep=""))]
mcmc.table.var <- as.data.frame(as.matrix(quartz(mcmc.var)))

quartz()
par(mfrow = c(1,2))
boxplot(mcmc.table.var[,1:18], at = 1:18, xlim=c(1,19), outline = FALSE, col = c(rep("white",times=9),"red","red",rep("white",times=7)))
}

# Boxplot N
# --------------------------------------------


# Will only work if dim theta and N > 1

if (length(as.vector(mcmc.table$'N[1]'))!=0){
  
var = "N"
mcmc.var = mcmc[,which(substr(varnames(mcmc),1,nchar(var))==paste(var,"",sep=""))]
mcmc.table.var <- as.data.frame(as.matrix(quartz(mcmc.var)))

boxplot(mcmc.table.var[,1:18], at = 1:18, xlim=c(1,19), outline = FALSE, col = c(rep("white",times=9),"red","red",rep("white",times=7)))

}


# ---------------------------------------------
# Joint and marginal posterior density
# ---------------------------------------------

# Only if parameters "a", "b" and  "c" exist in the model

if (length(as.vector(mcmc.table$'a'))!=0 & length(as.vector(mcmc.table$'b'))!=0 & length(as.vector(mcmc.table$'sigma'))!=0)
{
  
x1 <- as.vector(mcmc.table$'a')
x2 <- as.vector(mcmc.table$'b')
x3 <- as.vector(mcmc.table$'sigma')

def.par <- par(no.readonly = TRUE)
quartz()

par(pch='.')

pairs( cbind(x1,x2,x3),labels=c("a","b","sigma"), 
	 lower.panel=panel.smooth,
	 diag.panel=panel.dens, 
	 upper.panel=panel.cor,
	 cex.labels = 2.0, font.labels=1.2)
}

	 
# -----------------------------------------------
# Relationship theta ~ Q (or in the logit-scale)
# -----------------------------------------------

# Works only if "theta" and "logit_theta" is estimated for several years

if (length(as.vector(mcmc.table$'theta[1]'))!=0)
{
  
  var = "theta"
  mcmc.var = mcmc[,which(substr(varnames(mcmc),1,nchar(var))==paste(var,"",sep=""))]
  mean_theta <- summary(mcmc.var)$statistics[,1]
  
  quartz()
  par(mfrow = c(1,2))
  plot(data$Q,mean_theta[1:18], main = "theta as a function of Q", pch = 20, cex=1.5, col = c(rep("black",times=9),"red","red",rep("black",times=7)))
  
}

# Works only if "logit_theta" is estimated for several years

if (length(as.vector(mcmc.table$'logit_theta[1]'))!=0)
{
var = "logit_theta"
mcmc.var = mcmc[,which(substr(varnames(mcmc),1,nchar(var))==paste(var,"",sep=""))]
mean_logit_theta <- summary(mcmc.var)$statistics[,1]

plot(data$Q,mean_logit_theta[1:18], main = "logit(theta) as a function of Q", pch = 20, cex=1.5, col = c(rep("black",times=9),"red","red",rep("black",times=7)))
}


# ----------------------------------------------
# Trap efficiency  = f(Q)
# ----------------------------------------------

# Works only if "theta" is estimated for several years
if (length(as.vector(mcmc.table$'theta[1]'))!=0)
{

quartz()

par(mfrow = c(1,1), bty="n", mar=c(6,5,1,1))

size.labels = 1
size.text = 1.3
box.size1 = 1
box.size2 = 0.5
col1 <- "grey90"
col2 <- "grey55"
col1 = c(rep(col1, times=data$n)) ; col1[c(10,11)] <- "red"
col2 = c(rep(col2, times=data$n))

var = "theta"
mcmc.var = mcmc[,which(substr(varnames(mcmc),1,nchar(var))==paste(var,"",sep=""))]
mcmc.table.var <- as.data.frame(as.matrix(quartz(mcmc.var)))
X1 <- mcmc.table.var[,1:18]

x <- data$Q
# x.plot <- data$Q_pred
x.plot <- 1:8

labels <- ""

x.min = 0 ; x.max = max(x.plot)+1
y.min = 0 ; y.max = 1
y.label = "trap. eff"
x.label = "water discharge"
title = ""

boxplot(	X1[,1], 
		xaxt = "n", yaxt = "n", 
		xlim = c(x.min,x.max), ylim = c(y.min,y.max), 
		at = x[1], outpch = NA, boxwex=box.size1, col = col1[1]) 
for (i in 2:data$n)
{
boxplot(X1[,i], xaxt = "n", yaxt = "n", at = x[i], add = T, outpch = NA, boxwex=box.size1, col = col1[i])
}
axis(side =1, at=signif(seq(x.min,x.max,length.out=11),digits=2), labels = T, las=1, cex.axis=size.labels)
axis(side =2, cex.axis=size.labels)
mtext(x.label, line = 3, side = 1, cex = size.text)
mtext(y.label, line = 3, side = 2, cex = size.text)
mtext(title, line = 1.5, side = 3, adj=0, cex = size.text)

}

# ----------------------------------------------
# Trap efficiency  +  Q.pred and prediction
# ----------------------------------------------

# Works only if "theta_pred" (as a function of Q_pred defined in the data) exists in the model

if (length(as.vector(mcmc.table$'theta_pred[1]'))!=0)
{

quartz()

par(mfrow = c(1,1), bty="n", mar=c(6,5,1,1))

size.labels = 1
size.text = 1.3
box.size1 = 1
box.size2 = 0.5
col1 <- "grey90"
col2 <- "pink"
col1 = c(rep(col1, times=data$n)) ; col1[c(10,11)] <- "red"
col2 = c(rep(col2, times=data$n))

var = "theta"
mcmc.var = mcmc[,which(substr(varnames(mcmc),1,nchar(var))==paste(var,"",sep=""))]
mcmc.table.var <- as.data.frame(as.matrix(quartz(mcmc.var)))
# X1 = theta
X1 <- mcmc.table.var[,1:data$n]
# X2 = theta_pred
X2 <- mcmc.table.var[,(data$n+1):(data$n+data$n_pred)]

conf1l <- boxplot(X2)$stats[1,]
conf1u <- boxplot(X2)$stats[5,]
conf2l <- boxplot(X2)$stats[2,]
conf2u <- boxplot(X2)$stats[4,]

c <- median(mcmc.table$a)
d <- median(mcmc.table$b)
x <- data$Q
x.plot <- data$Q_pred
x.plot.inv <- sort(x.plot,decreasing=TRUE)
logit.theta.plot <- c*((x.plot-mean(data$Q))/sd(data$Q))+d
theta.plot <- exp(logit.theta.plot)/(1+exp(logit.theta.plot))

labels <- ""

x.min = 0 ; x.max = max(x.plot)+1
y.min = 0 ; y.max = 1
y.label = "trap. eff"
x.label = "water discharge"
title = ""

boxplot(	X1[,1], 
		xaxt = "n", yaxt = "n", 
		xlim = c(x.min,x.max), ylim = c(y.min,y.max), 
		at = x[1], outpch = NA, boxwex=box.size1, col = col1[1]) 

polygon(c(x.plot,x.plot.inv),c(conf1l,conf1u[length(data$Q_pred):1]),col="grey80",border=NA)
polygon(c(x.plot,x.plot.inv),c(conf2l,conf2u[length(data$Q_pred):1]),col="grey60",border=NA)

for (i in 1:data$n)
{
boxplot(X1[,i], xaxt = "n", yaxt = "n", at = x[i], add = T, outpch = NA, boxwex=box.size1, col = col1[i])
}

points(x.plot, theta.plot, type="l",lty=1, lwd=2, col="red")

axis(side =1, at=signif(seq(x.min,x.max,length.out=11),digits=2), labels = T, las=1, cex.axis=size.labels)
axis(side =2, cex.axis=size.labels)
mtext(x.label, line = 3, side = 1, cex = size.text)
mtext(y.label, line = 3, side = 2, cex = size.text)
mtext(title, line = 1.5, side = 3, adj=0, cex = size.text)

}


# Boxplot Juveniles --> Smolts survival
# --------------------------------------------


# Will only work if dim surv

if (length(as.vector(mcmc.table$'surv[1]'))!=0){
  
  var = "surv"
  mcmc.var = mcmc[,which(substr(varnames(mcmc),1,nchar(var))==paste(var,"",sep=""))]
  mcmc.table.var <- as.data.frame(as.matrix(quartz(mcmc.var)))
  
  quartz()
  boxplot(mcmc.table.var[,1:17], at = 1:17, xlim=c(1,18), outline = FALSE, col = c(rep("white",times=9),"red","red",rep("white",times=7)))
  
}

