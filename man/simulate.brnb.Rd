% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brnb.R
\name{simulate.brnb}
\alias{simulate.brnb}
\title{Simulate Responses}
\usage{
\method{simulate}{brnb}(object, nsim = 1, seed = NULL, ...)
}
\arguments{
\item{object}{an object representing a fitted model.}

\item{nsim}{number of response vectors to simulate. Defaults to 1.}

\item{seed}{an object specifying if and how the random number
generator should be initialized; see \code{\link[=set.seed]{set.seed()}} for details.}

\item{...}{extra arguments to be passed to methods. Not currently
used.}
}
\description{
Simulate one or more responses from the distribution corresponding
to a fitted model \code{\link[=brnb]{"brnb"}} object.
}
\examples{
# Example in Saha, K., & Paul, S. (2005). Bias-corrected maximum
# likelihood estimator of the negative binomial dispersion
# parameter.  Biometrics, 61, 179--185.
#
# Frequency distribution of red mites on apple leaves.
nomites <- 0:8
noleaves <- c(70, 38, 17, 10, 9, 3, 2, 1, 0)
fit_glmnb <- MASS::glm.nb(nomites~1,link="identity",weights = noleaves)
fit_brnb <- brnb(nomites ~ 1, link = "identity", transformation = "inverse",
                 type = "ML",weights = noleaves)
## Let us simulate 10 response vectors
sim_glmnb <- simulate(fit_glmnb, nsim = 10, seed = 123)
sim_brnb <-  simulate(fit_brnb, nsim = 10, seed = 123)
# The results from glm.nb and brnb with type = "ML" are
# exactly the same
all.equal(sim_glmnb, sim_brnb, check.attributes = FALSE)
}
