library("gnm")
library("aod")
data("barley", package = "gnm")
data("ships", package = "MASS")

qlmodel1 <- glm(y ~ site + variety, family = quasi(link = "logit", variance = "mu(1-mu)"), data = barley)
qlmodel2 <- glm(y ~ site + variety, family = wedderburn(link = "logit"), data = barley)

brmodel1 <- update(qlmodel1, method = "brquasi", type = "AS_empirical_mean")

brmodel2 <- update(qlmodel2, method = "brquasi", type = "AS_empirical_mean")


## Ships data, McCullagh and Nelder (1989, p205)
ships <- within(ships, {
    year <- factor(year)
    period <- factor(period)
})
ml_ships <- glm(incidents ~ offset(log(service)) + type + year + period, data = ships, family = quasipoisson(), subset = service != 0)
br_ships <- update(ml_ships, method = "brquasi", type = "AS_empirical_mean")

## Crabs from Agresti (2015) 8.1.3, p 271
Crabs <- read.table("http://users.stat.ufl.edu/~aa/glm/data/Crabs.dat", header = TRUE)

ql_crabs <- glm(y ~ weight, family = quasi(link = "log", variance = "mu"), data = Crabs)
br_crabs <- update(ql_crabs, method = "brquasi", type = "AS_empirical_mean", start = c(coef(ql_crabs), 1))


## Simulation - barley
data("barley", package = "gnm")
library("betareg")
eps <- 1e-05
barley_e <- within(barley, {
                   y <- ifelse(y == 0, eps, y)
                   y <- ifelse(y == 1, 1 - eps, y)})
gen_model <- betareg(y ~ site + variety, data = barley_e)
qlmodel2 <- glm(y ~ site + variety, family = quasi(link = "logit", variance = "mu(1-mu)"), data = barley_e)

truth <- coef(gen_model)
truth["(phi)"] <- 1/(1 + truth["(phi)"])

gen_model <- enrich(gen_model)

nsimu <- 1000
simu_data <- gen_model$auxiliary_functions$simulate(nsim= nsimu, seed = 123)

coefs_ql <- coefs_br <- coefs_ml <- grad_ql <- grad_br <- matrix(NA, length(truth), nsimu)
for (j in seq.int(nsimu)) {
    dat <- within(barley_e, { y = simu_data[, j]})
    ql <- try(update(qlmodel2, data = dat, method = "brquasi", type = "ML"))
    cql <- if (inherits(ql, "try-error")) rep(NA, length(truth)) else coef(ql, model = "full")
    br <- try(update(qlmodel2, data = dat, method = "brquasi", type = "AS_empirical_mean",
                     lambda = 1e-06, slowit = 0.5, maxit = 1000, trace = FALSE, only_beta = TRUE))
    cbr <- if (inherits(br, "try-error")) rep(NA, length(truth)) else coef(br, model = "full")  
    coefs_br[, j] <- cbr
    coefs_ql[, j] <- cql
    grad_ql[, j] <- if (inherits(ql, "try-error")) rep(NA, length(truth)) else ql$grad
    grad_br[, j] <- if (inherits(br, "try-error")) rep(NA, length(truth)) else br$grad
    if (j %% 5 == 0) cat(j, "\n")
}

rowMeans(coefs_ql, na.rm = TRUE) - truth
rowMeans(coefs_br, na.rm = TRUE) - truth


## Simulation - gasoline yield
library("betareg")
data("GasolineYield", package = "betareg")
eps <- 1e-05
gen_model <- betareg(yield ~ batch + temp, data = GasolineYield)
qlmodel2 <- glm(yield ~ batch + temp, family = wedderburn(link = "logit"), data = GasolineYield)

truth <- coef(gen_model)
truth["(phi)"] <- 1/(1 + truth["(phi)"])

gen_model <- enrich(gen_model)

nsimu <- 1000
simu_data <- gen_model$auxiliary_functions$simulate(nsim= nsimu, seed = 123)

coefs_ql <- coefs_br <- coefs_ml <- grad_ql <- grad_br <- matrix(NA, length(truth), nsimu)
for (j in seq.int(nsimu)) {
    dat <- within(GasolineYield, { yield = simu_data[, j]})
    ql <- try(update(qlmodel2, data = dat, method = "brquasi", type = "ML"))
    cql <- if (inherits(ql, "try-error")) rep(NA, length(truth)) else coef(ql, model = "full")
    br <- try(update(qlmodel2, data = dat, method = "brquasi", type = "AS_empirical_mean",
                     lambda = 1e-06, slowit = 0.5, maxit = 1000, trace = FALSE, start = cql))
    cbr <- if (inherits(br, "try-error")) rep(NA, length(truth)) else coef(br, model = "full")  
    coefs_br[, j] <- cbr
    coefs_ql[, j] <- cql
    grad_ql[, j] <- if (inherits(ql, "try-error")) rep(NA, length(truth)) else ql$grad
    grad_br[, j] <- if (inherits(br, "try-error")) rep(NA, length(truth)) else br$grad
    if (j %% 5 == 0) cat(j, "\n")
}

rowMeans(coefs_ql, na.rm = TRUE) - truth
rowMeans(coefs_br, na.rm = TRUE) - truth


## Negative binomial using the parameterization with var(Y) = mu_i * (1 + psi) / psi;
## phi = (1 + psi) / psi
## see McCullagh and Nelder (1989, Section 6.2.3)
library("parallel")

## overdispersion parameter is phi
psi2phi <- function(psi) (1 + psi)/psi
## psi in (0, infty) (see latent model)
phi2psi <- function(phi) 1/(phi - 1)

## This is in beta, log(phi) parameterization
loglik <- function(pars, y, x) {
    p <- ncol(x)
    beta <- pars[-(p + 1)]
    eta <- drop(x %*% beta)
    mu <- exp(eta)
    phi <- exp(pars[(p + 1)])
    psi <- phi2psi(phi)
    theta <- mu * psi    
    out <- sum(lgamma(theta + y) - lgamma(theta) + y * log(mu) + theta * log(theta) - (theta + y) * log(theta + mu))
    if (is.na(out) | is.infinite(out)) {
        -Inf
    }
    else {
        out
    }
}

## This is in beta, phi parameterization
simulate_nb <- function(pars, x) {
    n <- nrow(x)
    p <- ncol(x)
    beta <- pars[-(p + 1)]
    eta <- drop(x %*% beta)
    mu <- exp(eta)
    phi <- pars[(p + 1)]
    psi <- phi2psi(phi)
    theta <- mu * psi    
    rnegbin(n, mu, theta)
}

## This returns results in beta, phi parameterization
get_ml <- function(start, y, x, use_optim = FALSE) {
    p <- ncol(x)
    if (use_optim) {
        out <- optim(start, fn = loglik, x = x, y = y, method = "BFGS", control = list(fnscale = -1))
        obj <- -out$value
    }
    else {
        out <- nlminb(start, objective = function(pars) -loglik(pars, x = x, y = y))
        obj <- -out$objective
    }    
    coefs <- out$par
    gradient <- try(grad(function(pars) loglik(pars, y = y, x = x), x = coefs))
    if (inherits(gradient, "try-error")) {
        gradient <- rep(NA_real_, length(coefs))
    }  
    coefs[p + 1] <- exp(coefs[p + 1])
    gradient[p + 1] <- gradient[p + 1] / coefs[p + 1]
    list(coefficients = coefs, loglik = obj, gradient = gradient)
}


nobs <- 10
## Model matrix
set.seed(123)
x <- cbind(1, rbinom(nobs, 1, 0.5), rexp(nobs, 2))

repl <- 1
x <- do.call("rbind", replicate(repl, x, simplify = FALSE))

## parameters
beta <- c(2, 1, -1)
psi <- 1/5
truth <- c(beta, psi2phi(psi))

## simu size
nsimu <- 10000

set.seed(123)
simu_data <- replicate(nsimu, simulate_nb(truth, x))
## Tests for generator
## apply(simu_data, 1, mean) - exp(x %*% truth[1:3])
## apply(simu_data, 1, var)/apply(simu_data, 1, mean) - truth[4]


## coefs_ql <- coefs_qln <- coefs_br <- matrix(NA, 4, nsimu)
## converged_br <- numeric(nsimu)
## converged_ql <- numeric(nsimu)
## converged_qln <- numeric(nsimu)
## for (j in seq.int(nsimu)) {
##     y <- simu_data[, j]
##     ql <- glm(y ~ - 1 + x, family = quasi(link = "log", variance = "mu"), method = "brquasi",
##               type = "ML", slowit = 0.5, maxit = 1000, epsilon = 1e-04, only_beta = FALSE, disp_factor = "n-p")
##     qln <- glm(y ~ - 1 + x, family = quasi(link = "log", variance = "mu"), method = "brquasi",
##                type = "ML", slowit = 0.5, maxit = 1000, epsilon = 1e-04, only_beta = FALSE, disp_factor = "n")
##     coefs_ql[, j] <- c(coef(ql), summary(ql)$dispersion)
##     coefs_qln[, j] <- c(coef(qln, model = "full"))
##     converged_ql[j] <- ql$converged
##     converged_qln[j] <- qln$converged
##     ##
##     br <- glm(y ~ - 1 + x, family = quasi(link = "log", variance = "mu"), method = "brquasi",
##               type = "AS_empirical_mean", slowit = 0.5, maxit = 1000, only_beta = FALSE, disp_factor = "n")
##     converged_br[j] <- br$converged
##     coefs_br[, j] <- c(coef(br, model = "full"))
##     ##
##     if (j %% 1000 == 0) {
##         cat(paste("sample", j),
##             paste("br iter", br$iter),
##             paste("br max abs(grad)", max(abs(br$grad))), sep = " | ", "\n")
##     }
## }


res <- mclapply(seq.int(nsimu), function(j) {
    y <- simu_data[, j]
    ql <- glm(y ~ - 1 + x, family = quasi(link = "log", variance = "mu"), method = "brquasi",
              type = "ML", slowit = 0.5, maxit = 1000, epsilon = 1e-04, only_beta = FALSE, disp_factor = "n-p")
    qln <- glm(y ~ - 1 + x, family = quasi(link = "log", variance = "mu"), method = "brquasi",
               type = "ML", slowit = 0.5, maxit = 1000, epsilon = 1e-04, only_beta = FALSE, disp_factor = "n")
    br <- glm(y ~ - 1 + x, family = quasi(link = "log", variance = "mu"), method = "brquasi",
              type = "AS_empirical_mean", slowit = 0.5, maxit = 1000, only_beta = TRUE, disp_factor = "n-p")
        brn <- glm(y ~ - 1 + x, family = quasi(link = "log", variance = "mu"), method = "brquasi",
              type = "AS_empirical_mean", slowit = 0.5, maxit = 1000, only_beta = FALSE, disp_factor = "n")
    if (j %% 100 == 0) {
        cat(paste("sample", j),
            paste("br iter", br$iter),
            paste("br max abs(grad)", max(abs(br$grad))), sep = " | ", "\n")
    }
    list(coefs_br = coef(br, model = "full"),
         coefs_brn = coef(brn, model = "full"),
         coefs_ql = coef(ql, model = "full"),
         coefs_qln = coef(qln, model = "full"),
         grad_br = br$grad,
         grad_brn = brn$grad,
         grad_ql = ql$grad,
         grad_qln = qln$grad,
         converged_ql = ql$converged,
         converged_qln = qln$converged,
         converged_br = br$converged,
         converged_brn = brn$converged)   
}, mc.cores = 4)

cbind_list_element <- function(x, name) do.call("cbind", lapply(x, "[[", name))
c_list_element <- function(x, name) do.call("c", lapply(x, "[[", name))

coefs_br <- cbind_list_element(res, "coefs_br")
coefs_brn <- cbind_list_element(res, "coefs_brn")
coefs_ql <- cbind_list_element(res, "coefs_ql")
coefs_qln <- cbind_list_element(res, "coefs_qln")

grad_br <- cbind_list_element(res, "grad_br")
grad_brn <- cbind_list_element(res, "grad_brn")
grad_ql <- cbind_list_element(res, "grad_ql")
grad_qln <- cbind_list_element(res, "grad_qln")

converged_br <- c_list_element(res, "converged_br")
converged_brn <- c_list_element(res, "converged_brn")
converged_ql <- c_list_element(res, "converged_ql")
converged_qln <- c_list_element(res, "converged_qln")

converged <- converged_br & converged_ql & converged_qln & converged_brn

rowMeans(coefs_ql[, converged]) - truth
rowMeans(coefs_qln[, converged]) - truth
rowMeans(coefs_br[, converged]) - truth
rowMeans(coefs_brn[, converged]) - truth

rowMeans((coefs_ql - truth)^2)
rowMeans((coefs_qln - truth)^2)
rowMeans((coefs_br - truth)^2)
rowMeans((coefs_brn - truth)^2)

sqrt(rowMeans((coefs_ql - truth)^2)/ nsimu)
