#' Auxiliary function for \code{\link{glm}} fitting using the
#' \code{\link{brglmFit}} method.
#'
#' Typically only used internally by \code{\link{brglmFit}}, but may
#' be used to construct a \code{control} argument.
#'
#' @inheritParams stats::glm.control
#' @param epsilon positive convergence tolerance epsilon
#' @param maxit integer giving the maximal number of iterations
#'     allowed
#' @param trace logical indicating if output should be produced for
#'     each iteration
#' @param type the type of fitting methodo to be used. The options are
#'     \code{AS-mean} (mean-bias reducing adjusted scores; default),
#'     \code{correction} (asymptotic bias correction) and
#'     \code{ML} (maximum likelihood).
#' @param transformation the transformation of the dispersion to be
#'     estimated. Default is \code{identity}. See Details.
#' @param slowit a positive real used as a multiplier for the
#'     stepsize. The smaller it is the smaller the steps are
#' @param max_step_factor the maximum number of step halving steps to
#'     consider
#'
#' @details \code{\link{brglmControl}} provides default values and
#'     sanity checking for the various constants that control the
#'     iteration and generally the behaviour of
#'     \code{\link{brglmFit}}.
#'
#'      When \code{trace} is true, calls to \code{cat} produce the
#'      output for each iteration.  Hence, \code{options(digits = *)}
#'      can be used to increase the precision.
#'
#'      \code{transformation} sets the transformation of the
#'      dispersion parameter for which the bias reduced estimates are
#'      computed. Can be one of "identity", "sqrt", "inverse", "log"
#'      and "inverseSqrt". Custom transformations are accommodated by
#'      supplying a list of two expressions (transformation and
#'      inverse transformation). See the examples for more details.
#'
#'
#' @return a list with components named as the arguments, including
#'     symbolic expressions for the dispersion transformation
#'     (\code{Trans}) and its inverse (\code{inverseTrans})
#'
#' @author Ioannis Kosmidis \email{i.kosmidis@ucl.ac.uk}
#'
#' @seealso \code{\link{brglmFit}} and \code{\link{glm.fit}}
#'
#' @examples
#'
#' data("coalition", package = "brglm2")
#' ## The maximum likelihood fit with log link
#' coalitionML <- glm(duration ~ fract + numst2, family = Gamma, data = coalition)
#'
#' ## Bias reduced estimation of the dispersion parameter
#' coalitionBRi <- update(coalitionML, method = "brglmFit")
#' coef(coalitionBRi, model = "dispersion")
#'
#' ## Bias reduced estimation of log(dispersion)
#' coalitionBRl <- update(coalitionML, method = "brglmFit",  transformation = "log")
#' coef(coalitionBRl, model = "dispersion")
#'
#' ## Just for illustration: Bias reduced estimation of dispersion^0.25
#' my_transformation <- list(expression(dispersion^0.25), expression(transformed_dispersion^4))
#' coalitionBRc <- update(coalitionML, method = "brglmFit",
#'                        transformation = my_transformation)
#' coef(coalitionBRc, model = "dispersion")
#'
brglmControl <- function(epsilon = 1e-10, maxit = 100,
                         trace = FALSE,
                         type = c("AS_mean", "correction", "ML"),
                         transformation = "identity",
                         slowit = 1,
                         max_step_factor = 12) {
    type <- match.arg(type)

    if (is.character(transformation)) {
        Trans <- switch(transformation,
                        identity = expression(dispersion),
                        sqrt = expression(dispersion^0.5),
                        inverse = expression(1/dispersion),
                        log = expression(log(dispersion)),
                        inverseSqrt = expression(1/sqrt(dispersion)),
                        stop(transformation, " is not one of the implemented dispersion transformations"))
        inverseTrans <- switch(transformation,
                               identity = expression(transformed_dispersion),
                               sqrt = expression(transformed_dispersion^2),
                               inverse = expression(1/transformed_dispersion),
                               log = expression(exp(transformed_dispersion)),
                               inverseSqrt = expression(1/transformed_dispersion^2))
    }
    else {
        if (is.list(transformation) && (length(transformation) == 2)) {
            Trans <- transformation[[1]]
            inverseTrans <- transformation[[2]]
            transformation <- "custom_transformation"
        }
        else {
            stop("transformation can be either one of 'identity', 'sqrt', 'inverse', 'log' and 'inverseSqrt', or a list of two expressions")
        }
    }
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of 'epsilon' must be > 0")
    list(epsilon = epsilon, maxit = maxit, trace = trace,
         type = type,
         Trans = Trans,
         inverseTrans = inverseTrans,
         transformation = transformation,
         slowit = slowit,
         max_step_factor = max_step_factor)
}

