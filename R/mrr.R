#' Calculate median rate ratio (MRR) from a variance
#'
#' mor: Calculate median odds ratio (MOR) from a random effect variance from
#' a binomial GLMM. inv.mor: Inverse function to get variance from MOR.
#' For a poisson GLMM these functions will transform variances to median
#' *rate* ratios (MRR), and vice versa. mrr and inv.mrr are aliases
#' for mor and inv.mor.
#' See:
#'   Interpreting Parameters in the Logistic Regression Model with Random Effects
#'   Author(s): Klaus Larsen, J¯rgen Holm Petersen, Esben Budtz-J¯rgensen, Lars Endahl
#'   Biometrics, Vol. 56, No. 3 (Sep., 2000), pp. 909-914
#' Equations:
#'      mor=exp(sqrt(2*v)*qnorm(0.75))   (MOR function)
#' =>   log(mor)=sqrt(2*v)*qnorm(0.75)
#' =>   (log(mor)/qnorm(0.75))^2=2*v
#' =>   v=((log(mor)/qnorm(0.75))^2)/2   (inverse MOR function)
#'
#' @param v Variance, input to mor and mrr
#' @export
#' @examples
#' # a random effect variance of 1.3 between levels (e.g. sites)...
#' mor(1.3)
#' # ...implies that a typical pair of randomly chosen levels
#' # will differ in odds (or rate, for count GLMMs) by a factor
#' # of 3.
#' # inverse function
#' inv.mor(3)
#' # mrr and inv.mrr are aliases for mor and inv.mor.
#' mrr(1.3)
#' inv.mrr(3)

mrr<-function(v)exp(sqrt(2*v)*qnorm(0.75))
