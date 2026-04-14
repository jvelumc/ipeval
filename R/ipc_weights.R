flip_surv_event <- function(formula) {
  lhs <- formula[[2]]

  lhs[[3]] <- call("!", lhs[[3]])   # negate event indicator
  formula[[2]] <- lhs
  formula
}



#' Shift censored observation times to break ties with events
#'
#' Adjusts observed times for censored observations by adding a small epsilon,
#' ensuring that censored cases are treated as occurring slightly after events
#' when tied on the recorded time scale.
#'
#' @details In survival data, event times are often recorded with limited
#'   precision, leading to ties between failure events and censoring times.
#'
#'   A common convention in Cox partial likelihood computations is that, at a
#'   tied time point, events are assumed to occur before censoring. Many
#'   standard packages for estimating Cox models are built around this
#'   assumption.
#'
#'   When estimating the censoring distribution (rather than the event
#'   distribution) with a Cox model, the usual approach is to model
#'   \code{Surv(time, status == 0) ~ ... }. However, this flips the definition
#'   of censored cases and events, and now the censored observations (which are
#'   handled as events) are treated to occur before the events.
#'
#'   With the \code{prodlim} package the censoring distribution can be computed
#'   with \code{prodlim::prodlim(..., reverse = TRUE)}, correctly handling ties.
#'   An implementation for Cox does not exist as far as I'm aware. As a crude
#'   workaround, this function can be used to shift censored observation times
#'   by some small amount such that there are no more ties left. A way to
#'   compute the censoring distribution with Cox while correctly taking to
#'   account ties between censored cases and events would be better.
#'
#' @param time Numeric vector of follow-up times.
#' @param status Integer or logical event indicator (1 = event, 0 = censored).
#' @param epsilon Small positive perturbation added to censored times to break
#'   ties.
#'
#' @returns A numeric vector of modified times with censored observations
#'   shifted forward by \code{epsilon}.
#' @export
#'
#' @examples
#' time <- c(5, 5, 8, 10)
#' status <- c(1, 0, 1, 0)
#' push_censored_observations(time, status)
push_censored_observations <- function(time, status, epsilon = 1e-8) {
  time[status == 0] <- time[status == 0] + epsilon
  return(time)
}


ipc_weights <- function(data, formula, type, time_horizon) {
  # perhaps use riskRegression::ipcw()
  if (type == "KM")
    stopifnot(rhs_is_one(formula))

  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)

  p_uncensored <- switch(
    type,
    KM = {
      fit <- prodlim::prodlim(formula, data = data, reverse = TRUE)
      p_not_censor <- stats::stepfun(fit$time, c(1, fit$surv), right = TRUE)
      list(
        model = fit,
        probability = p_not_censor(pmin(y[, "time"], time_horizon))
      )
    },
    cox = {
      # coxph has no reverse argument, need to flip it manually
      # does not handle event/censor ties correctly
      # throw a warning if there are ties
      flipped_form <- flip_surv_event(formula)

      fit <- survival::coxph(flipped_form, data = data, model = TRUE, x = TRUE)

      event_times <- unique(y[y[, "status"] == 1, "time"])
      censor_times <- unique(y[y[, "status"] == 0, "time"])
      if (any(event_times %in% censor_times)) {
        warning("There are events that have the same time as censored
                observations. Consider running ?push_censored_observations().")
      }

      list(
        model = fit,
        probability = 1-predict_cox(fit, data, pmin(y[, "time"], time_horizon))
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
