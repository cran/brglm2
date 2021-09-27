library("brglm2")
load("inst/franscisco/franscisco.Rdata")

fm <- Mig ~ State + Gender + SC + SPI + FL + SQ + DM +
    DP + HS + EP + WSS + HBE + ACLD + BC + Transfer + SS + Agegrp +
    Gender * Agegrp

mod <- glm(fm, family = "binomial", data = df_All,
               method = "brglm_fit",
               slowit = 0.5)

mod0 <- glm(fm, family = "binomial", data = df_All,
            method = "brglm_fit", trace = 1)

summary(mod)
