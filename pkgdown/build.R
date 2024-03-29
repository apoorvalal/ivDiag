
setwd("~/github/ivDiag")

# initializing
library(usethis)
library(sinew)
library(pkgdown)
usethis::use_readme_rmd()
usethis::use_pkgdown()
usethis::use_news_md() # update logs

# remember to knitr README.Rmd
devtools::install_github("r-lib/pkgdown")
pkgdown::build_site(install = FALSE)

# or alternatively
setwd("~/github/ivDiag")
library(pkgdown)
#init_site()
build_home()
build_news()
build_reference()
build_articles()
