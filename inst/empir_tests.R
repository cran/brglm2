library("numDeriv")
lizardsBR_mean <- glm(cbind(grahami, opalinus) ~ height + diameter +
                      light + time, family = binomial(probit), data = lizards,
                      method = "dispersion")

## brglmFit
theta <- theta
all.equal(unname(umat(theta, dispersion_block = TRUE)),
          unname(numDeriv:::hessian(function(pars) gradient(pars, only_beta = FALSE)[4], x = theta)),
          tolerance = 1e-06)

## beta
uu <- umat(theta)
sapply(1:3, function(jj) {
    all.equal(unname(numDeriv::hessian(function(pars) gradient(pars, only_beta = TRUE)[jj], x = theta)),
              unname(uu[[jj]]), tolerance = 1e-10)
})


## j mat

## beta
sapply(1:3, function(j) {
    all.equal(unname(grad(function(pars) gradient(pars)[j], x = theta)),
              unname(-jmat(theta)[j, ]), tolerance = 1e-10)
})

## dispersion
all.equal(unname(grad(function(pars) gradient(pars, only_beta = FALSE)[4], x = theta)),
          unname(-jmat(theta, only_beta = FALSE)[4, ]), tolerance = 1e-10)

## dmat
par_settings <- expand.grid(r = 1:4, s = 1:4, t = 1:4)

out1 <- apply(par_settings, 1, function(x) {
    r <- x[1]
    s <- x[2]
    t <- x[3]
    dpsii <- sapply(1:nobs, function(i) grad(function(pars) gradient(pars, only_beta = FALSE, contributions = TRUE)[i, r], x = theta))[s, ]
    psii <- gradient(theta, only_beta = FALSE, contributions = TRUE)[, t]
    cat(r, s, t, "\n")
    sum(dpsii * psii)
})


out1 <- as.list(numeric(4))
for (r in 1:4) {
    dm <- matrix(NA, 4, 4)
    for (s in 1:4) {
        for (t in 1:4) {
            dpsii <- sapply(1:nobs, function(i) grad(function(pars) gradient(pars, only_beta = FALSE, contributions = TRUE)[i, r], x = theta))[s, ]
            psii <- gradient(theta, only_beta = FALSE, contributions = TRUE)[, t]
            cat(r, s, t, "\n")
            dm[s, t] <- sum(dpsii * psii)
        }
    }
    out1[[r]] <- dm
}

mats <- c(dmat(theta, only_beta = FALSE), list(dmat(theta, dispersion_block = TRUE, only_beta = FALSE)))
out2 <- apply(par_settings, 1, function(x) {
    r <- x[1]
    s <- x[2]
    t <- x[3]
    oo <- mats[[r]][s, t]
    names(oo) <- paste0("r", r, "s", s, "t", t)
    oo
})
all.equal(out1, out2, tolerance = 1e-08)

##
qmu <- function(pars) {
    fit <- key_quantities(pars, y = y, level = 2, qr = FALSE)
    with(fit, {
        emus <- (y - mus)
        bmus <- weights * d1mus / varmus
        d1bmus <- weights * (d2mus / varmus - d1mus^2 * d1varmus / varmus^2)
        d2bmus <- weights * (d3mus / varmus - d2mus * d1mus * d1varmus / varmus^2 - 2 * d1mus * d2mus * d1varmus / varmus^2 + 2 * d1mus^3 * d1varmus^2 / varmus^3)
        d2bmus <- weights * (d3mus / varmus - (3 * d1mus * d2mus * d1varmus + d1mus^3 * d2varmus) / varmus^2 + 2 * d1mus^3 * d1varmus^2 / varmus^3)
        list(bmus = bmus,
             d1bmus = d1bmus,
             d2bmus = d2bmus)
        ## list(a = d1mus^2 * d1varmus / varmus^2,
        ##      b = ((2 * d1mus * d2mus * d1varmus + d1mus^3 * d2varmus) * varmus^2 - 2 * d1mus^3 * d1varmus^2 * varmus)/varmus^4)
        
        ## list(qmus = (bmus * d1mus - d1bmus * emus) / dispersion,
        ##      d1qmus = (2 * d1bmus * d1mus + bmus * d2mus - d2bmus * emus) / dispersion)
    })
}
