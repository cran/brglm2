## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6
)

## ---- echo = TRUE-------------------------------------------------------------
library("brglm2")
data("stemcell", package = "brglm2")

## ---- echo = TRUE-------------------------------------------------------------
stemcells_ml <- bracl(research ~ as.numeric(religion) + gender, weights = frequency, data = stemcell,
                      parallel = TRUE, type = "ML")
summary(stemcells_ml)

## ---- echo = TRUE-------------------------------------------------------------
class(stemcells_ml)

## ---- echo = TRUE-------------------------------------------------------------
stemcells_ml_full <- bracl(research ~ as.numeric(religion) + gender, weights = frequency, data = stemcell,
                           parallel = FALSE, type = "ML")
summary(stemcells_ml_full)

## ---- echo = TRUE-------------------------------------------------------------
(lrt <- deviance(stemcells_ml) - deviance(stemcells_ml_full))

## ---- echo = TRUE-------------------------------------------------------------
(df1 <- df.residual(stemcells_ml) - df.residual(stemcells_ml_full))

## ---- echo = TRUE-------------------------------------------------------------
pchisq(lrt, df1, lower.tail = FALSE)

## ---- echo = TRUE-------------------------------------------------------------
summary(update(stemcells_ml, type = "AS_mean"))

## ---- echo = TRUE-------------------------------------------------------------
summary(update(stemcells_ml, type = "AS_median"))

## ---- echo = TRUE-------------------------------------------------------------
predict(stemcells_ml, type = "probs")

