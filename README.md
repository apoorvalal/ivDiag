# ivDiag

R package for IV diagnostics accompanying Lal, Lockhart, Xu, and Zu (2021). Computes bootstrap SEs, F-stats, effective F stats, and Local-to-zero IV estimates.

## Installation

```r
library(remotes)
install_github("apoorvalal/ivDiag")
```

---

### Notes (to self/coauthors) for site maintenance 

Examples are all currently in vignette. If we want examples for each
documentation entry,

+ Add/edit examples in source files in `R/*.R`
+ Run `./setup.R` to rebuild documentation
  + (optionally) edit `_pkgdown.yml` and `DESCRPTION` for aesthetics + metadata
+ Run `./build_site.R` to rebuild site
