pkgname <- "ggcustomalluvial"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ggcustomalluvial')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("normalize_alluvium_sizes")
### * normalize_alluvium_sizes

flush(stderr()); flush(stdout())

### Name: normalize_alluvium_sizes
### Title: Normalize a Vector to Proportions
### Aliases: normalize_alluvium_sizes

### ** Examples

library(dplyr)
df <- tibble(group = c("A", "B", "C"), size = c(3, 6, 1))



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
