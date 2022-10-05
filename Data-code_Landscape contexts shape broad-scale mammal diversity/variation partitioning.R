#################### variation partitioning #####################

# take lm_1km_local and lm_1km_context10km as an example

# input species richness data (dependent Variable)
sr <- read.csv("species richness.csv")
SR <- sr$SR
# input explanatory variables
local_raw <- read.csv("lm_1km_local.csv") # local landscape
context_raw <- read.csv("lm_1km_context10km.csv") # landscape context
macro_raw <- read.csv("macro-environmental variables.csv") # macro-environment

#------------------------------------------------------------------------------------------
##### variable selection
# prescreen local landscape, landscape context and macro-environment variables based on VIF (variance inflation factor), respectively

# take local landscape variables as an example
library(car)
library(MASS)
Id <- local_raw$Id # ID
local_raw <- local_raw[!names(local_raw) %in% c("Id")] 
# remove perfect multicolliniarity (remove explanatory variables with coefficient = NA)
ols <- lm(SR ~., data = local_raw)
coef <- summary(ols)[["coefficients"]]
coefs <- coef[-1,1]
coefs <- data.frame(gsub('[local_raw]', '', names(coefs)), coefs)
colnames(coefs) <- c("col", "coefficients")
local_raw <- local_raw[,colnames(local_raw) %in% coefs$col]
# calculate variance inflation factor (VIF) of variables
ols <- function(local_raw){
  local_raw <- as.data.frame(scale(local_raw))
  ols <- lm(SR ~ ., data = local_raw)
  vifs <- data.frame(names(vif(ols)), vif(ols))
  colnames(vifs) <- c("col", "vif")
  vifs
}
vifs <- ols(local_raw)
# remove explanatory variables with VIF > 5 iteratively
while(max(vifs$vif) > 5) {
  vifs <- vifs[-which(vifs$vif==max(vifs$vif)), ]
  local_raw <- local_raw[,colnames(local_raw) %in% vifs$col]
  vifs <- ols(local_raw)
}
local <- local_raw[,colnames(local_raw) %in% vifs$col] # local landscape variables after selection
local <- cbind(Id, local, SR) # add Id & SR



#------------------------------------------------------------------------------------------
############ variation partitioning

# combine explanatory variables (after variable selection)
lc <- cbind(local, context)
cm <- cbind(context, macro)
lm <- cbind(local, macro)
lcm <- cbind(local, context, macro)

# define the function to calculate adjusted-R2 based on the ordinary least squares (OLS) model
rsqfunc <- function(var) {
    SR <- var$SR
    # remove the column of Id and SR to retain only explanatory variables
    var <- var[, !names(var) %in% c("Id", "SR")]
    var <- as.matrix(scale(var))
    # specify OLS model
    model <- lm(SR ~ var)
    # output adjusted-R2
    rsq <- summary(model)$adj.r.squared
    rsq
  }
   
# define the function to calculate pseudo-R2 based on the random forests (RF) model
library(randomForest)
rsqfunc <- function(var){
    SR <- var$SR
    # remove the column of Id and SR to retain only explanatory variables
    var <- var[,!names(var)%in%c("Id", "SR")]
    var <- as.matrix(scale(var))
    # specify RF model (by default, mtry = 1/3 of total number of variables)
    RF <- randomForest(SR ~ ., data = var, ntree = 500)
    # output pseudo-R2, defined as 1-[mean square error/var(y)]
    rsq <- RF$rsq[500]
    rsq
  }

    
### calculate adjusted-R2 (based on the OLS model) or pseudo-R2 (based on the RF model) with bootstrapping

sampsize <- ceiling(nrow(sr)/2) # random sample size = half of total number of grid cells
nrep <- 999 # number of replications for bootstrapping
set.seed(123)

rsqAll_local <- matrix(nrow = nrep, ncol = 1)
rsqAll_context <- matrix(nrow = nrep, ncol = 1)
rsqAll_macro <- matrix(nrow = nrep, ncol = 1)
rsqAll_lc <- matrix(nrow = nrep, ncol = 1)
rsqAll_cm <- matrix(nrow = nrep, ncol = 1)
rsqAll_lm <- matrix(nrow = nrep, ncol = 1)
rsqAll_lcm <- matrix(nrow = nrep, ncol = 1)

for (i in 1:nrep) {
  # random sampling  
  randsp <- sample(1:nrow(local), sampsize, replace = FALSE)
  rspdata_local <- local[randsp, ]
  rspdata_context <- context[randsp, ]
  rspdata_macro <- macro[randsp, ]
  rspdata_lc <- lc[randsp, ]
  rspdata_cm <- cm[randsp, ]
  rspdata_lm <- lm[randsp, ]
  rspdata_lcm <- lcm[randsp, ]
  # calculate adjusted-R2 (based on the OLS model) or pseudo-R2 (based on the RF model)
  rsq_local <- rsqfunc(rspdata_local)
  rsq_context <- rsqfunc(rspdata_context)
  rsq_macro <- rsqfunc(rspdata_macro)
  rsq_lc <- rsqfunc(rspdata_lc)
  rsq_cm <- rsqfunc(rspdata_cm)
  rsq_lm <- rsqfunc(rspdata_lm)
  rsq_lcm <- rsqfunc(rspdata_lcm)
  # output adjusted-R2 (based on the OLS model) or pseudo-R2 (based on the RF model)
  rsqAll_local[i, ] <- rsq_local
  rsqAll_context[i, ] <- rsq_context
  rsqAll_macro[i, ] <- rsq_macro
  rsqAll_lc[i, ] <- rsq_lc
  rsqAll_cm[i, ] <- rsq_cm
  rsqAll_lm[i, ] <- rsq_lm
  rsqAll_lcm[i, ] <- rsq_lcm
}


##### Independent explanatory power of landscape context and local landscape

L <- rsqAll_lcm - rsqAll_cm # independent explanatory power of local landscape
C <- rsqAll_lcm - rsqAll_lm # independent explanatory power of landscape context
# calculate mean and standard deviation (in percentage)
L_mn <- mean(L) * 100 
C_mn <- mean(C) * 100
L_sd <- sd(L) * 100
C_sd <- sd(C) * 100


##### calculate variation independently and shared explained by local landscape, landscape context and macro-environment

L <- rsqAll_lcm - rsqAll_cm # variation independently explained by local landscape
C <- rsqAll_lcm - rsqAll_lm # variation independently explained by landscape context
M <- rsqAll_lcm - rsqAll_lc # variation independently explained by macro-environment
LC <- rsqAll_lm + rsqAll_cm - rsqAll_macro - rsqAll_lcm # variation shared explained by local landscape and landscape context
CM <- rsqAll_lc + rsqAll_lm - rsqAll_local - rsqAll_lcm # variation shared explained by landscape context and macro-environment
LM <- rsqAll_lc + rsqAll_cm - rsqAll_context - rsqAll_lcm # variation shared explained by local landscape and macro-environment
LCM <- rsqAll_local + rsqAll_context + rsqAll_macro - rsqAll_lc - rsqAll_cm - rsqAll_lm + rsqAll_lcm # variation shared explained by local landscape, landscape context and macro-environment   
vpt <- cbind(L, C, M, LC, CM, LM, LCM)
# calculate mean and standard deviation of variation explained (in percentage)
ve_mn <- apply(vpt, 2, mean) * 100
ve_sd <- apply(vpt, 2, sd) * 100
