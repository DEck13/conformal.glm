pkgname <- "paraconformal"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('paraconformal')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("conformalprediction")
### * conformalprediction

flush(stderr()); flush(stdout())

### Name: conformalprediction
### Title: Prediction Regions for Generalized Linear Regression Models
### Aliases: conformalprediction
### Keywords: generalized linear regression models conformal prediction

### ** Examples

# example of section 2.4 in Geyer (2009)
# data(sports)
# out <- glmdr(cbind(wins, losses) ~ 0 + ., family = "binomial", data = sports)
#summary(out)



cleanEx()
nameEx("regions")
### * regions

flush(stderr()); flush(stdout())

### Name: regions
### Title: Prediction Regions for Generalized Linear Regression Models
### Aliases: regions
### Keywords: generalized linear regression models conformal prediction

### ** Examples

# example of section 2.4 in Geyer (2009)
# data(sports)
# out <- glmdr(cbind(wins, losses) ~ 0 + ., family = "binomial", data = sports)
#summary(out)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
