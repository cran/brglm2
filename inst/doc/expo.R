## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6
)

## ----echo = TRUE--------------------------------------------------------------
library("brglm2")
data("aids", package = "brglm2")
aids

## ----echo = TRUE--------------------------------------------------------------
aids_mod <- glm(cbind(symptomatic, asymptomatic) ~ AZT + race, 
                  family = binomial(), data = aids)
summary(aids_mod)

## ----echo = TRUE--------------------------------------------------------------
expo(aids_mod, type = "ML")

## ----echo = TRUE--------------------------------------------------------------
expo(aids_mod, type = "correction*")
expo(aids_mod, type = "Lylesetal2012")
expo(aids_mod, type = "correction+")

## ----echo = TRUE--------------------------------------------------------------
expo(aids_mod, type = "AS_median")

## ----echo = TRUE--------------------------------------------------------------
data("endometrial", package = "brglm2")
endometrialML <- glm(HG ~ NV + PI + EH, data = endometrial, family = binomial())
endometrialML

## ----echo = TRUE--------------------------------------------------------------
library("detectseparation")
update(endometrialML, method = detect_separation)

## ----echo = TRUE--------------------------------------------------------------
expo(endometrialML, type = "correction*")
expo(endometrialML, type = "correction+")
expo(endometrialML, type = "Lylesetal2012")

## ----echo = TRUE--------------------------------------------------------------
aids_mod_br <- update(aids_mod, method = "brglmFit")
expo(aids_mod_br, type = "correction*")

