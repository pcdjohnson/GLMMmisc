#' Check fit of (g)lmer residuals-v-fitted plot by overlaying simulated points
#'
#' Function to assess the fit of a GLMM by making a residuals-v-fitted-values plot
#' and overlaying residuals and fitted values from from a model fitted to data
#' simulated from the fitted model. The rationale is that, although we often don't
#' know how a resid-v-fitted plot should look for a GLMM, we do know that if we
#' simulate from the fitted model, then refit the original model to the simulated data,
#' the resulting fit will be perfect. If the fit is perfect, then the residuals and fitted
#' values from the simulated fit should also be perfect (but with sampling error).
#' Therefore, if we make resid-v-fitted plots from both the real and simulated fitted
#' models, we can judge the fit of the model on how similar they look. Works on lme4 model
#' objects, not tested on others yet.
#'
#' @param object A model object fitted with glmer or lmer
#' @param add.sim.resid Overlay simulated residuals? If FALSE, only the real residuals are plotted.
#' @export
#' @examples
#' # fit a Poisson-lognormal GLMM (a Poisson GLMM with an
#' # observation-level random effect [OLRE])
#' # to the grouseticks data (see ?grouseticks)
#' fit.poisln <-
#'   glmer(TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD) + (1 | LOCATION) + (1 | INDEX),
#'         family = "poisson", data = grouseticks)
#'
#' # plot the Pearson residuals v fitted values
#' plotloess(residfitted.olre(fit.poisln), loess.col = "black")
#'
#' # when we assess a residuals-v-fitted-values plot we are looking for
#' # there to be no trend (in particular no nonlinear trend) and
#' # homoscedasticity (homogeneous variance of the residuals along the
#' # fitted values -- no fanning or tapering).
#' # this isn't always easy to assess by eye, e.g. here the spread of the points
#' # seems to get narrower at higher fitted values, but is that a sign of decreasing
#' # variance or just decreasing density?
#'
#' # often people look for a normal distribution of Pearson residuals as well...
#' hist(residfitted.olre(fit.poisln)$PearsonResiduals) # not normal
#' # ...but in fact there is no expectation of normality unless the error distribution
#' # itself is normal.
#'
#' # an aside... when the model is Poisson-lognormal (i.e. has an OLRE), simply using
#' # plot(fit.poisln) gives the wrong residuals and fitted values
#' # ("wrong" in the sense that they generally won't be trendless and homoscedastic even if the
#' # model fits perfectly):
#' plot(fit.poisln)
#' # ...for reasons explained here
#' # https://github.com/pcdjohnson/miscR/blob/master/residfitted.olre.R
#'
#' # returning to the original plot
#' plotloess(residfitted.olre(fit.poisln), loess.col = "black")
#'
#' # is that how the plot should look, if the model fits well?
#' # with GLMMs, especially non-gaussian ones, it's often very hard to tell,
#' # because how the plot should look is often dependent on characteristics of
#' # the model and the data.
#'
#' # compare real and simulated resid vs fitted scatter plots.
#' sim.residplot(fit.poisln)
#' # (repeat this line a few times as there will be variation between
#' # simulated data sets.)
#'
#' # the Poisson-lognormal seems to fit pretty well.
#' # try fitting a Gaussian GLMM to the tick counts.
#' # just by looking at the distribution of the tick counts we would
#' # guess that a model with Gaussian errors will fit badly
#' hist(grouseticks$TICKS)
#'
#' # fit the model (minus the overdispersion random effect)
#' fit.gauss <-
#'   lmer(TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD) + (1 | LOCATION),
#'        data = grouseticks)
#'
#' # plot residuals X fitted values
#' plot(fit.gauss)
#'
#' # ...doesn't look good.
#' # in the case of a Gaussian GLMM we not only expect the plot to show
#' # a lack of trend and homoscedasticy, we also expect the residuals to
#' # be normally distributed around zero. clearly these last two assumptions
#' # are violated here.
#'
#' # try comparing with a simulated (and perfectly fitting) model
#' sim.residplot(fit.gauss, add.sim.resid = FALSE)
#' sim.residplot(fit.gauss)
#' # remember to repeat the plot a few times to get an idea of the influence
#' # of sampling error.
#'
#' # the plot clearly illustrates how badly fitting the Gaussian model is.
#'
#' # try a binomial GLMM, where the response is binary: >0 ticks counted.
#' grouseticks$TICKS.any <- grouseticks$TICKS > 0.1
#' table(grouseticks$TICKS.any)
#'
#' fit.binom <-
#'   glmer(TICKS.any ~ YEAR + scale(HEIGHT) + (1 | BROOD) + (1 | LOCATION),
#'         family = "binomial", data = grouseticks)
#'
#' # plot residuals X fitted values
#' plot(fit.binom)
#'
#' # this sort of weird pattern is typical of binomial GL(M)Ms fitted to binary data.
#' # it's very hard to detect absence or presence of a trend or homoscedasticity.
#'
#' # try comparing with a simulated (and perfectly fitting) model
#' sim.residplot(fit.binom)
#' # that looks very similar.
#'
#' # another binomial example, taken from ?glmer
#' cbpp$obs <- 1:nrow(cbpp)
#' (gm2 <- glmer(cbind(incidence, size - incidence) ~ period +
#'                 (1 | herd) +  (1|obs),
#'               family = binomial, data = cbpp))
#' plot(gm2)
#' # this is the kind of ugly trend we typically get when a binomial GLMM has an
#' # observation-level random effect (the "obs" effect here).
#' # (I haven't yet worked out how to fix this within the residfitted.olre
#' # function, but see
#' # https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q1/021818.html).
#' # the advantage of simulating from the model is that we don't have to get rid
#' # of this trend. if the trend is a sign of poor fit, it won't appear in
#' # the simulated scatter plot...
#' sim.residplot(gm2) # repeat a few times
#' # ...but it does, so suggesting that the residuals aren't showing signs of poor fit.



sim.residplot <-
  function(object, add.sim.resid = TRUE) {
    require(lme4)
    if(!add.sim.resid) {
      plotloess(residfitted.olre(object, warn = FALSE),
                 loess.col = 1)
    } else {
      object.sim <- refit(object, simulate(object)[[1]])
      xy <-
        cbind(residfitted.olre(object, warn = FALSE),
              residfitted.olre(object.sim, warn = FALSE))
      names(xy) <- paste0(names(xy), c("", "", ".sim", ".sim"))
      plotloess(residfitted.olre(object, warn = FALSE),
                 xlim = range(xy$Fitted, xy$Fitted.sim),
                 ylim = range(xy$PearsonResiduals, xy$PearsonResiduals.sim),
                 loess.col = 1)
      xy$loess.line.sim <- predict(loess(PearsonResiduals.sim ~ Fitted.sim, data = xy))
      xy <- xy[order(xy$Fitted.sim), ]
      points(xy[,c("Fitted.sim", "PearsonResiduals.sim")], pch = 3, col = 2)
      lines(xy[,c("Fitted.sim", "loess.line.sim")], type = "l", col = 2)
      y.all <- unlist(xy[, grep("Residuals", names(xy))])
      x.all <- unlist(xy[, grep("Fitted", names(xy))])
      y <- y.all[y.all > (max(y.all) * 2/3)]
      x <- x.all[y.all > (max(y.all) * 2/3)]
      leg.breaks <- c(-Inf, min(x.all) + diff(range(x.all)) * c(1, 2)/3, Inf)
      leg.pos <- c("topleft", "top", "topright")[which.min(table(cut(x, leg.breaks)))]
      legend(leg.pos, legend = c("Real", "Simulated"), pch = c(1, 3), col = 1:2, bty = "n")
    }
  }


