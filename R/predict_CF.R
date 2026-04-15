predict_CF <- function(model, data, A_column, CF_treatment, time_horizon) {
  # predict outcome probabilities for all patients, setting their treatment
  # to CF_treatment
  data[[A_column]] <- CF_treatment
  if ("glm" %in% class(model)) {
    return(predict_glm(model, data))
  }
  if ("coxph" %in% class(model)) {
    return(predict_cox(model, data, time_horizon))
  }
  stop("model class ", class(model), " not supported")
}


predict_glm <- function(model, data) {
  unname(stats::predict(model, newdata = data, type = "response"))
}


predict_cox <- function(model, data, time_horizon) {

  bh <- survival::basehaz(model, centered = FALSE)

  n <- nrow(data)

  if (length(time_horizon) == 1L) {
    time_horizon <- rep(time_horizon, n)
  } else if (length(time_horizon) != n) {
    stop("time_horizon must have length 1 or nrow(data)")
  }

  idx <- findInterval(time_horizon, bh$time, left.open = TRUE)

  H0_horizon <- numeric(n)
  valid <- idx > 0
  H0_horizon[valid] <- bh$hazard[idx[valid]]

  S0_horizon <- exp(-H0_horizon)

  lp <- stats::predict(model, newdata = data, type = "lp")

  S_horizon <- S0_horizon^exp(lp)
  1-S_horizon

}


