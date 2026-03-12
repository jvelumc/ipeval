ipc_weights <- function(data, formula, type, time_horizon) {


  if (type == "KM")
    stopifnot(rhs_is_one(formula))

  mf <- model.frame(formula, data)
  y <- model.response(mf)



  p_uncensored <- switch(
    type,
    KM = {
      fit <- prodlim::prodlim(formula, data = data, reverse = TRUE)
      p_not_censor <- stepfun(fit$time, c(1, fit$surv), right = TRUE)
      list(
        model = fit,
        probability = p_not_censor(pmin(y[, "time"], time_horizon))
      )
    },
    cox = {
      # coxph has no reverse argument, need to flip it manually
      if ("cfscore_time" %in% colnames(mf) || "cfscore_status" %in% colnames(mf)) {
        stop("Please don't use the variable names cfscore_time or cfscore_status in your data")
      }
      mf$cfscore_time <- y[, "time"]
      mf$cfscore_status <- y[, "status"]

      flipped_form <- update.formula(formula, Surv(cfscore_time, cfscore_status == 0) ~ .)

      fit <- coxph(flipped_form, data = mf, model = TRUE, x = TRUE)
      list(
        model = fit,
        probability = 1-predict_cox(fit, mf, pmin(time_horizon, mf$cfscore_time))
      )
    },
    stop("cens.model ", type, " not implemented")
  )

  # if censored before time horizon, weight is 0,
  # else, weight is 1/probability uncensored at event/time horizon
  w <- ifelse(
    y[, "status"] == 0 & y[, "time"] < time_horizon,
    0,
    1 / p_uncensored$probability
  )

  list(
    model = p_uncensored$model,
    weights = unname(w)
  )
}
