bootstrap_iteration <- function(data, ip_score) {
  # this could probably be refactored a bit, but will require efforts to
  # refactor the main ipscore function too

  # works by creating a new ipscore object based on the original, where all
  # required data has been resampled (& new ipt & ipc weights)

  bs_sample <- sample(nrow(data), size = nrow(data), replace = T)
  bs_ipt <- ipt_weights(data[bs_sample, ], ip_score$ipt$propensity_formula)
  bs_pred <- lapply(ip_score$predictions, function(x) x[bs_sample])

  bs_ipscore <- list(
    outcome_type = ip_score$outcome_type,
    metrics = ip_score$metrics,
    predictions = bs_pred,
    outcome = ip_score$outcome[bs_sample],
    observed_treatment = ip_score$observed_treatment[bs_sample],
    treatment_of_interest = ip_score$treatment_of_interest
  )
  bs_ipscore$ipt$weights <- bs_ipt$weights

  # are we using stabilized weights?
  if (!is.null(ip_score$ipt$stable_model)) {
    stable_treatment_formula <- stats::update.formula(ip_score$ipt$propensity_formula, . ~ 1)
    bs_sipt <- ipt_weights(data[bs_sample, ], stable_treatment_formula)
    bs_ipscore$ipt$weights <- 1/bs_sipt$weights * bs_ipt$weights
  }


  if (ip_score$outcome_type == "survival") {
    bs_ipc <- ipc_weights(
      data = data[bs_sample, ],
      formula = ip_score$ipc$cens_formula,
      type = ip_score$ipc$method,
      time_horizon = ip_score$time_horizon
    )
    bs_ipscore$status_at_horizon <- ip_score$status_at_horizon[bs_sample]
    bs_ipscore$ipc$weights <- bs_ipc$weights
  }
  metrics <- get_metrics(bs_ipscore)

  return(metrics)
}


bootstrap <- function(data, ip_score) {
  b <- lapply_progress(
    as.list(1:ip_score$bootstrap_iterations),
    function(x) {
      bootstrap_iteration(data, ip_score)
    },
    "bootstrapping",
    progress = ip_score$bootstrap_progress
  )
  # transpose results
  # (iteration > metric > model) -> (metric > model > iteration)

  # for calibration plot:
  # (iteration > metric > [pred/obs, model]) ->
  # (metric > model > iteration > list(pred = , obs = ))
  transposed <- lapply(ip_score$metrics, function(m) {
    P <- lapply(names(ip_score$predictions), function(p) {
      if (m != "calplot") { # 1 numeric result, simple to combine & transpose
        sapply(b, function(i) i[[m]][[p]])
      } else { # calibration plot, consisting of 2 vectors of preds & obs
        lapply(b, function(i) {
          list(
            pred = i[[m]][["pred", p]],
            obs = i[[m]][["obs", p]]
          )
        })
      }
    })
    names(P) <- names(ip_score$predictions)
    P
  })
  names(transposed) <- ip_score$metrics

  # # summarize
  conf.int <- lapply(ip_score$metrics, function(m) {
    CI <- lapply(names(ip_score$predictions), function(p) {
      if (m != "calplot") {
        return(ci(transposed[[m]][[p]], cover = 0.95))
      } else {
        return(NA)
      }
    })
    names(CI) <- names(ip_score$predictions)
    CI
  })
  names(conf.int) <- ip_score$metrics

  list(
    results = conf.int,
    raw = transposed
  )
}


lapply_progress <- function(x, FUN, task_description, progress = TRUE) {
  # same as lapply, but print a progress indicator
  n <- length(x)

  if (progress == FALSE) {
    result <- lapply(as.list(1:n), FUN)
  } else {
    FUN2 <- function(x, i, n) {
      result <- FUN(x)
      cat("\r", task_description, ": ", i, "/", n, "     ")
      return(result)
    }

    result <- lapply(as.list(1:n), function(i) FUN2(x[[i]], i, n))
    cat("\r")
  }

  return(result)
}

ci <- function(values, cover = 0.95) {
  lower <- (1-cover) / 2
  upper <- 1 - lower
  stats::quantile(values, probs = c(lower, upper))
}
