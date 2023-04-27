setwd('/Users/xyq/GitHub/ivDiag/vignettes')
load('/Users/xyq/GitHub/ivDiag/vignettes/GSZ2016.rds')

head(gsz_2016)
gsz <- as.data.frame(gsz_2016)
haven::zap_labels(gsz)
class(gsz)

library(ivDiag)
data(ivDiag)

Y <- "totassoc_p"
D <- "libero_comune_allnord"
Z <- "bishopcity"
weights <- "population"
cl <- NULL
controls <- c('altitudine', 'escursione', 'costal', 'nearsea', 'population', 
              'pop2', 'gini_land', 'gini_income')

gsz <- gsz[, unique(c(Y, D, Z, controls, FE, cl, weights))]
gsz <- as.data.frame(gsz[complete.cases(gsz), ])
dim(gsz)

data <- as.data.frame(readRDS("rueda2017.rds"))
Y <- "e_vote_buying" # Y: outcome of interest
D <-"lm_pob_mesa" # D: endogenous treatment
Z <- "lz_pob_mesa_f" # Z: instrumental variable
controls <- c("lpopulation", "lpotencial") # covariates of control variables
cl <- "muni_code" # clusters
weights <- FE <- NULL # no weights or fixed effects
data <- data[, unique(c(Y, D, Z, controls, FE, cl, weights))]
data <- haven::zap_labels(data)
rueda <- as.data.frame(data[complete.cases(data), ])
dim(rueda)

save(rueda, gsz, gsz_south, file = "ivDiag.RData")

print(g <- ivDiag(data = gsz, Y = Y, D = D, Z = Z, controls = controls, weights = weights))

zfs_prior # estimate in ZFS south sample - saved in processing

gsz_south <- haven::zap_labels(haven::read_dta("gsz2016_south.dta"))
gsz_south <- gsz_south[, unique(c(Y, Z, controls, FE, cl, weights,"capoluogo"))]
gsz_south <- as.data.frame(gsz_south[complete.cases(gsz_south), ])
dim(gsz_south)


# 0.178, 0.137




# vanilla clustered fit

m1 = felm(f, data = gsz_2016, weights = gsz_2016$population) 

source("/Users/xyq/GitHub/ivDiag/R/ltz.R")
source("/Users/xyq/GitHub/ivDiag/R/utils.R")
ltz_est <- ltz(gsz_2016, Y, D, Z, controls, weights = weights, zfs.coef = 0.177, zfs.se = 0.137)
ltz_est







library(ivDiag);
library(magrittr);
library(patchwork);
library(ggfortify);
library(lfe)

source("/Users/xyq/GitHub/ivDiag/R/viz_iv_dists.R")



