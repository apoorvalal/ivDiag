library(ivDiag); library(data.table); library(lfe)

# %% bootstrap iv test
n = 1e3
z = rnorm(n)
d = rbinom(n, 1, plogis(z))
y = 2 * d + rnorm(n)
df = data.frame(y, d, z)
# %% bootstrap
boot_IV(df, 'y', 'd', 'z', nboots = 100)

# %%
dgp_zfs = function(γ = 0, n = 1e4, share_zfs = 0.2,
    α_0 = -1, π_0 = 2, β_0 = 1, σ_β = 0.25,  sed = 42){
  Z = rbinom(n, 1, 0.5)
  α = rnorm(n, α_0, 1)
  π = runif(n, π_0 - .5, π_0 + .5)
  e = rnorm(n, 0, 1)
  β = rnorm(n, β_0, σ_β)
  d = data.table(Z, α, π, β, e)
  # zero first stage population - last (1-s) obs
  d[round((1 - share_zfs)*n):n, zfs_pop := 1][zfs_pop == 1, π := 0]
  # generate data
  d[,
    D_star := α + π*Z + e][,
    D := ifelse(D_star > 0, 1, 0)][,
    Y := γ*Z + β*D + e]
}
# %%
d = dgp_zfs(γ = 0.2, β_0 = 2)
# zero first stage - nonzero
zfs_fs = felm(D ~ Z, d[!is.na(zfs_pop)])
# reduced form
zfs_rf = felm(Y ~ Z, d[!is.na(zfs_pop)])
# %%
tsls = felm(Y ~ 1 | 0 | (D ~ Z), d)
# %%
ltz(d, 'D', 'Z', zfs_fs$coefficients[2], 0.05, tsls)
