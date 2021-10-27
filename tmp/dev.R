# %% ####################################################
rm(list = ls())
library(LalRUtils)
LalRUtils::libreq(tidyverse, data.table, fst, fixest, rio, foreach,
                  janitor, tictoc, RColorBrewer, patchwork, RPushbullet, IRdisplay)
theme_set(lal_plot_theme()) # add _d() for dark
set.seed(42)
# %% load all scripts
library(ivDiag)
# lapply(list.files("../R", full.names = TRUE), source)
# ls()
# %%
df <-readRDS("apsr_Meredith_2013.rds")
Y <- "DemShareDB"; D <- "DemShareGOV"; Z <- "HomeGOV"
controls <- c("HomeDB")
cl <- "fips"
FE <- c("fips","RaceID")
weights <- "Weight"

# %%
bootres=boot_IV(data=df, Y=Y, D=D, Z=Z, controls=controls, FE = FE,
  cl=cl, weights=weights, nboot = 50, sens = T)
bootres |> print()
# %%
