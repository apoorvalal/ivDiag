# ivDiag

R package for IV diagnostics accompanying Lal, Lockhart, Xu, and Zu (2021). Computes bootstrap SEs, F-stats, effective F stats, and Local-to-zero IV estimates.

## Installation

```r
library(remotes)
install_github("apoorvalal/ivDiag")
```

### Bootstrap

```r
n = 1e3
z = rnorm(n); ε = rnorm(n)
d = rbinom(n, 1, plogis(z + ε))
y = 2 * d + ε
df = data.frame(y, d, z)
# biased
lm(y ~ d)
# %% bootstrap
boot_IV(df, 'y', 'd', 'z', nboots = 100)

```

```
Bootstrapping:
Parallelising 100 reps on 7 cores
Bootstrap took 7.079 sec.

$est_ols
Coef                2.6846  SE.t0.0571  SE.b0.0525  CI.b 2.5%2.583  CI.b 97.5%2.775
$est_2sls
Coef                2.016   SE.t0.1533  SE.b0.1561  CI.b 2.5%1.67   CI.b 97.5%2.298
$F_stat
F.standard    187.0336      F.robust    241.9001    F.cluster<NA>   F.boot   205.8122
$p_iv 1
$N 1000
$N_cl NULL
```

### Zero-first stage test

CHR 2012 'Plausibly Exogenous' estimates with prior derived from known
`zero-first-stage` subsample.

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
