## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6
)

## ----echo = TRUE--------------------------------------------------------------
freq <- c(15, 16, 16, 27, 33, 20,
          21, 18, 26, 41, 38, 27,
          29, 21, 33, 60, 41, 42)
dose <- rep(c(0, 10, 33, 100, 333, 1000), 3)
plate <- rep(1:3, each = 6)
(salmonella <- data.frame(freq, dose, plate))

## ----echo = TRUE--------------------------------------------------------------
ames_f <- freq ~ dose + log(dose + 10)

## ----echo = TRUE--------------------------------------------------------------
library("brglm2")
ames_ML <- brnb(ames_f, link = "log", data = salmonella,
                transformation = "identity",  type = "ML")
## Estimated regression and dispersion parameters
est <- coef(ames_ML, model = "full")
## Estimated standard errors for the regression parameters
sds <- sqrt(c(diag(ames_ML$vcov.mean), ames_ML$vcov.dispersion))
round(cbind(est, sds), 4)

## ----echo = TRUE--------------------------------------------------------------
ames_BC <- update(ames_ML, type = "correction")
## Estimated regression and dispersion parameters
est <- coef(ames_BC, model = "full")
## Estimated standard errors for the regression parameters
sds <- sqrt(c(diag(ames_BC$vcov.mean), ames_BC$vcov.dispersion))
round(cbind(est, sds), 4)

## ----echo = TRUE--------------------------------------------------------------
ames_BRmean <- update(ames_ML, type = "AS_mean")
## Estimated regression and dispersion parameters
est <- coef(ames_BRmean, model = "full")
## Estimated standard errors for the regression parameters
sds <- sqrt(c(diag(ames_BRmean$vcov.mean), ames_BRmean$vcov.dispersion))
round(cbind(est, sds), 4)

## ----echo = TRUE--------------------------------------------------------------
ames_BRmedian <- update(ames_ML, type = "AS_median")
## Estimated regression and dispersion parameters
est <- coef(ames_BRmedian, model = "full")
## Estimated standard errors for the regression parameters
sds <- sqrt(c(diag(ames_BRmedian$vcov.mean), ames_BRmedian$vcov.dispersion))
round(cbind(est, sds), 4)

## ----echo = TRUE--------------------------------------------------------------
ames_BRmixed <- update(ames_ML, type = "AS_mixed")
## Estimated regression and dispersion parameters
est <- coef(ames_BRmixed, model = "full")
## Estimated standard errors for the regression parameters
sds <- sqrt(c(diag(ames_BRmixed$vcov.mean), ames_BRmixed$vcov.dispersion))
round(cbind(est, sds), 4)

