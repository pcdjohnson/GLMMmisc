#' Output "corrected" Pearson residuals and fitted values
#'
#' Output "corrected" Pearson residuals and fitted values from a Poisson
#' GLMM with an observation-level random effect (OLRE) fitted using lme4::glmer.
#' Defaults to residuals(., type = "pearson") and fitted(.) for other classes of fit.
#'
#' @param object Fitted model object
#' @param warn Warn when not Poisson glmer object with OLRE
#' @export
#' @examples
#' library(lme4)
#' form1 <- TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD) + (1 | LOCATION) + (1 | INDEX)
#' fit.poisln  <- glmer(form1, family = "poisson", data = grouseticks)
#' residfitted.olre(fit.poisln)

residfitted.olre <-
  function(object, warn = TRUE) {
    require(AICcmodavg)
    require(lme4)
    response<-model.frame(object)[[1]]
    n <- length(response)
    re <- lme4::ranef(object)
    re.length <- sapply(re, nrow)
    od.term <- names(re)[re.length == n]
    pois.log.norm <-
      length(od.term) == 1 && fam.link.mer(object)$family == "poisson" && fam.link.mer(object)$link == "log"
    if(!pois.log.norm) {
      if(warn) warning("Fitted model is not Poisson-lognormal. Using stats::residuals and stats::fitted.")
      f <- fitted(object)
      r <- residuals(object, type = "pearson")
    }
    if(pois.log.norm) {
      od.ranef <- re[[od.term]][[1]]
      f <- exp(log(fitted(object)) - od.ranef)
      r <- (response - f) / sqrt(f + (f^2) * c(exp(lme4::VarCorr(object)[[od.term]]) - 1))
    }
    return(data.frame(Fitted = f, PearsonResiduals = r))
  }
