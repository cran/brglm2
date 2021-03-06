#' brglm2: Bias Reduction in Generalized Linear Models
#'
#' Estimation and inference from generalized linear models using
#' implicit and explicit bias reduction methods (Kosmidis, 2014), and
#' other penalized maximum likelihood methods. Currently supported
#' methods include the mean bias-reducing adjusted scores approach in
#' Firth (1993) and Kosmidis & Firth (2009), the median bias-reduction
#' adjusted scores approach in Kenne Pagui et al. (2017), the
#' correction of the asymptotic bias in Cordeiro & McCullagh (1991),
#' the mixed bias-reduction adjusted scores approach in Kosmidis et al
#' (2020), maximum penalized likelihood with powers of the Jeffreys
#' prior as penalty, and maximum likelihood.
#'
#'
#' In the special case of generalized linear models for binomial,
#' Poisson and multinomial responses (both nominal and ordinal), mean
#' and median bias reduction and maximum penalized likelihood return
#' estimates with improved frequentist properties, that are also
#' always finite, even in cases where the maximum likelihood estimates
#' are infinite (e.g. complete and quasi-complete separation in
#' multinomial regression; see also \code{\link{detect_separation}}
#' and \code{\link{check_infinite_estimates}} for pre-fit and post-fit
#' methods for the detection of infinite estimates in binomial
#' response generalized linear models). Estimation in all cases takes
#' place via a modified Fisher scoring algorithm, and S3 methods for
#' the construction of confidence intervals for the reduced-bias
#' estimates are provided.
#'
#' The core model fitters are implemented by the functions
#' \code{\link{brglm_fit}} (univariate generalized linear models),
#' \code{\link{brmultinom}} (baseline category logit models for
#' nominal multinomial responses), and \code{\link{bracl}} (adjacent
#' category logit models for ordinal multinomial responses).
#'
#' @details
#'
#'
#' The similarly named **brglm** R package can only handle generalized
#' linear models with binomial responses. Special care has been taken
#' when developing **brglm2** in order not to have conflicts when the
#' user loads **brglm2** and **brglm** simultaneously. The development
#' and maintenance of the two packages will continue in parallel,
#' until **brglm2** incorporates all **brglm** functionality and gets
#' an appropriate wrapper to the \code{brglm::brglm} function.
#'
#' @author Ioannis Kosmidis \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso
#'
#' \code{\link{brglm_fit}}, \code{\link{brmultinom}}, \code{\link{bracl}}
#'
#' @references
#'
#' Kosmidis I, Firth D (2020). Jeffreys-prior penalty, finiteness
#' and shrinkage in binomial-response generalized linear
#' models. *Biometrika* \doi{10.1093/biomet/asaa052}
#'
#' Cordeiro G M, McCullagh P (1991). Bias correction in generalized
#' linear models. *Journal of the Royal Statistical Society. Series B
#' (Methodological)*, **53**, 629-643 \doi{10.1111/j.2517-6161.1991.tb01852.x}
#'
#' Firth D (1993). Bias reduction of maximum likelihood estimates,
#' Biometrika, **80**, 27-38 \doi{10.2307/2336755}
#'
#' Kenne Pagui E C, Salvan A, Sartori N (2017). Median bias
#' reduction of maximum likelihood estimates. *Biometrika*, **104**,
#' 923–938 \doi{10.1093/biomet/asx046}
#'
#' Kosmidis I, Kenne Pagui E C, Sartori N (2020). Mean and median bias
#' reduction in generalized linear models. *Statistics and Computing*,
#' **30**, 43-59 \doi{10.1007/s11222-019-09860-6}
#'
#' Kosmidis I, Firth D (2009). Bias reduction in exponential family
#' nonlinear models. *Biometrika*, **96**, 793-804 \doi{10.1093/biomet/asp055}
#'
#' Kosmidis I, Firth D (2010). A generic algorithm for reducing
#' bias in parametric estimation. *Electronic Journal of Statistics*,
#' **4**, 1097-1112 \doi{10.1214/10-EJS579}
#'
#' Kosmidis I (2014). Bias in parametric estimation: reduction and
#' useful side-effects. *WIRE Computational Statistics*, **6**,
#' 185-196 \doi{10.1002/wics.1296}
#'
#' @docType package
#' @name brglm2
#' @import stats
#' @import enrichwith
#' @import Matrix
#' @import MASS
#' @importFrom graphics plot
#' @importFrom nnet class.ind
#' @importFrom numDeriv grad
#'
#'
NULL

## Suggestion by Kurt Hornik to avoid a warning related to the binding
## of n which is evaluated by family$initialize
if (getRversion() >= "2.15.1") globalVariables(c("n", "lambda"))

#' Generic method for checking for infinite estimates
#' @param object a fitted model object (e.g. the result of a
#'     \code{\link{glm}} call).
#' @param ... other options to be passed to the method.
#'
#'
#' @note
#'
#'
#' \code{check_infinite_estimates} will be removed from \pkg{brglm2}
#' at version 0.8. An new version of
#' \code{check_infinite_estimates} is now maintained in the
#' \pkg{detectseparation} R package at
#' \url{https://cran.r-project.org/package=detectseparation}.
#'
#' @seealso check_infinite_estimates.glm
#'
#' @export
check_infinite_estimates <- function(object, ...) {
    function_moves_to_new_package(gsub("\\(|\\)", "", deparse(match.call()[1])),
                                  "0.8",
                                  "brglm2",
                                  "detectseparation")
    UseMethod("check_infinite_estimates")
}

