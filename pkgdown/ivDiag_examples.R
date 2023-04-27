### Sample Code for *ivDiag* ###


##########################
## Installation
##########################

library(remotes)
install_github("apoorvalal/ivDiag")


##########################
## Example 1: Rueda (2017)
##########################

rm(list=ls())
library(ivDiag)
data(ivDiag)
ls()

Y <- "e_vote_buying" # Y: outcome of interest
D <-"lm_pob_mesa" # D: endogenous treatment
Z <- "lz_pob_mesa_f" # Z: instrumental variable
controls <- c("lpopulation", "lpotencial") # covariates of control variables
cl <- "muni_code" # clusters
weights <- FE <- NULL # no weights or fixed effects

# first stage (raw)
par(mar = c(4, 4, 2, 2))
plot(rueda$lz_pob_mesa_f, rueda$lm_pob_mesa, col = "#777777", cex = 0.5, 
     main = "Raw Data", xlab = "Instrument", ylab = "Treatment")
abline(lm(lm_pob_mesa ~ lz_pob_mesa_f, data = rueda), col = 2, lwd = 2, lty = 2)

# first stage (partial out)
z_res <- lm(lz_pob_mesa_f ~ lpopulation + lpotencial, data = rueda)$residuals
d_res <- lm(lm_pob_mesa ~ lpopulation + lpotencial, data = rueda)$residuals
plot(z_res, d_res, col = "#777777", cex = 0.5, main = "Covariates Partialled Out", 
     xlab = "Residualized Instrument", ylab = "Residualized Treatment")
abline(lm(d_res ~ z_res), col = 2, lwd = 2, lty = 2)

# omnibus funcdtion
library(ivDiag)
g <- ivDiag(data = rueda, Y=Y, D = D, Z = Z, controls = controls, cl = cl)
names(g)
g

# plot coefficients and CIs
plot_coef(g)


## seperate functions
eff_F(data = rueda, Y = Y, D = D, Z = Z, controls = controls, cl = cl, 
      FE = NULL, weights = NULL)

AR_test(data = rueda, Y = Y, D = D, Z = Z, controls = controls, cl = cl, 
        FE = NULL, weights = NULL)

tF(coef = -0.9835, se = 0.1540, Fstat = 8598)

##########################
## Example 2: GSZ (2016)
##########################

Y <- "totassoc_p"
D <- "libero_comune_allnord"
Z <- "bishopcity"
weights <- "population"
controls <- c('altitudine', 'escursione', 'costal', 'nearsea', 'population', 'pop2', 'gini_land', 'gini_income')

# table instrument and treatment
table(gsz$bishopcity)
table(gsz$libero_comune_allnord)
table(gsz$bishopcity, gsz$libero_comune_allnord)

# distribution of the outcome
hist(gsz$totassoc_p, breaks = 100, xlab = "#Nonprofit Organizations Per Capita", main = "")


## Omnibus function
g <- ivDiag(data = gsz, Y = Y, D = D, Z = Z, controls = controls, weights = weights)
g
plot_coef(g)


## Zero-first-stage
library(estimatr)
zfs <- lm_robust(totassoc_p ~ bishopcity + altitudine + escursione + costal + nearsea + 
          population + pop2 + gini_land + gini_income + capoluogo, data = gsz_south, 
          weights = gsz_south$population, se_type = "HC1")

summary(zfs)$coefficients["bishopcity", 1:2]

## Local-to-zero adjustment
ltz_est <- ltz(data = gsz, Y = Y, D = D, Z = Z, controls = controls, weights = weights, 
              mu = 0.178, Sig = 0.137^2)
ltz_est

## visualization
viz_iv_dists(g, ltz_est, xlim = c(-0.5, 7)) 

