# ivDiag

R package for IV diagnostics accompanying [Lal, Lockhart, Xu, and Zu (2021)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3905329). Computes bootstrap SEs, F-stats, effective F stats, and Local-to-zero IV estimates.

## Suggested Citation

Lal, Apoorva and Lockhart, Mackenzie William and Xu, Yiqing and Zu, Ziwen, How Much Should We Trust Instrumental Variable Estimates in Political Science? Practical Advice based on Over 60 Replicated Studies (August 14, 2021). Available at SSRN: https://ssrn.com/abstract=3905329 or http://dx.doi.org/10.2139/ssrn.3905329

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
