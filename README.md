# ivDiag (in progress)

R package for IV diagnostics accompanying Lal, Lockhart, Xu, and Zu (2021). Computes bootstrap SEs, F-stats, effective F stats, and Local-to-zero IV estimates. 


## Bootstrap

```r
n = 1e3
z = rnorm(n)
d = rbinom(n, 1, plogis(z))
y = 2 * d + rnorm(n)
df = data.frame(y, d, z)
# %% bootstrap
boot_IV(df, 'y', 'd', 'z', nboots = 100)
```

```
Bootstrapping:
Parallelising 100 reps on 7 cores 
Bootstrap took 8.245 sec.
$est_ols
Coef1.9645SE.t0.0646SE.b0.0706CI.b 2.5%1.834CI.b 97.5%2.107
$est_2sls
Coef2.2018SE.t0.152SE.b0.1544CI.b 2.5%1.876CI.b 97.5%2.441
$F_stat
F.standard223.9963F.robust329.0763F.cluster<NA>F.boot311.826
$p_iv
1
$N
1000
$N_cl
NULL
```

## Zero-first stage test (CHR 'Plausibly Exogenous' test)

```r
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
# zero first stage
zfs_fs = felm(D ~ Z, d[!is.na(zfs_pop)])
zfs_fs |> summary()
# reduced form
zfs_rf = felm(Y ~ Z, d[!is.na(zfs_pop)])
zfs_rf |> summary()
# %%
tsls = felm(Y ~ 1 | 0 | (D ~ Z), d)
# %%
ltz(d, 'D', 'Z', zfs_fs$coefficients[2], 0.05, tsls)

```

```
$est
2.522
$se
0.344
$ci
1.847614628891753.19611942854956
```
