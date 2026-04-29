acl_probs <- function(theta, X, ncat, po = TRUE) {
    X <- as.matrix(X)
    int <- (ncat - 1):1
    if (po) {
        nslopes <- ncol(X) - 1L
        if (length(theta) != (ncat - 1) + nslopes) {
            stop("for `po = TRUE`, `length(theta)` must be `ncat - 1 + ncol(X) - 1`")
        }
        coefs <- cbind(rev(cumsum(theta[int])),
                       int * matrix(theta[-int], nrow = ncat - 1,
                                    ncol = nslopes, byrow = TRUE))
    } else {
        if (length(theta) != (ncat - 1) * ncol(X)) {
            stop("for `po = FALSE`, `length(theta)` must be `(ncat - 1) * ncol(X)`")
        }
        coefs <- matrix(theta, nrow = ncat - 1)
        coefs <- apply(coefs, 2, function(x) rev(cumsum(x[int])))
        if (is.null(dim(coefs))) {
            coefs <- matrix(coefs, ncol = 1L)
        }
    }
    eta <- X %*% t(coefs)
    fits <- cbind(eta, 0)
    exp_fits <- exp(fits - apply(fits, 1, max))
    exp_fits / rowSums(exp_fits)
}


loglik <- function(theta, Y, X, po = TRUE) {
    X <- as.matrix(X)
    if (is.factor(Y)) {
        y <- as.integer(Y)
        ncat <- nlevels(Y)
    } else {
        ylev <- sort(unique(Y))
        y <- match(Y, ylev)
        ncat <- length(ylev)
    }
    probs <- acl_probs(theta, X, ncat = ncat, po = po)
    sum(log(pmax(probs[cbind(seq_along(y), y)], .Machine$double.xmin)))
}
