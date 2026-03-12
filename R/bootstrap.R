bootstrap_iteration <- function(data, cfscore) {
  # this could probably be refactored a bit, but will require efforts to
  # refactor the main cfscore function too

  # works by creating a new cfscore object based on the original, where all
  # required data has been resampled (& new ipt & ipc weights)

  bs_sample <- sample(nrow(data), size = nrow(data), replace = T)
  bs_ipt <- ipt_weights(data[bs_sample, ], cfscore$ipt$propensity_formula)
  bs_pred <- lapply(cfscore$predictions, function(x) x[bs_sample])

  bs_cfscore <- list(
    outcome_type = cfscore$outcome_type,
    metrics = cfscore$metrics,
    predictions = bs_pred,
    outcome = cfscore$outcome[bs_sample],
    observed_treatment = cfscore$observed_treatment[bs_sample],
    treatment_of_interest = cfscore$treatment_of_interest
  )
  bs_cfscore$ipt$weights <- bs_ipt$weights

  # are we using stabilized weights?
  if (!is.null(cfscore$ipt$stable_model)) {
    stable_treatment_formula <- update.formula(cfscore$ipt$propensity_formula, . ~ 1)
    bs_sipt <- ipt_weights(data[bs_sample, ], stable_treatment_formula)
    bs_cfscore$ipt$weights <- 1/bs_sipt$weights * bs_ipt$weights
  }


  if (cfscore$outcome_type == "survival") {
    bs_ipc <- ipc_weights(
      data = data[bs_sample, ],
      formula = cfscore$ipc$cens_formula,
      type = cfscore$ipc$method,
      time_horizon = cfscore$time_horizon
    )
    bs_cfscore$status_at_horizon <- cfscore$status_at_horizon[bs_sample]
    bs_cfscore$ipc$weights <- bs_ipc$weights
  }
  metrics <- get_metrics(bs_cfscore)

  return(metrics)
}


bootstrap <- function(data, cfscore) {
  b <- lapply_progress(
    as.list(1:cfscore$bootstrap_iterations),
    function(x) {
      bootstrap_iteration(data, cfscore)
    },
    "bootstrapping",
    progress = cfscore$bootstrap_progress
  )
  # transpose results
  # (iteration > metric > model) -> (metric > model > iteration)

  # for calibration plot:
  # (iteration > metric > [pred/obs, model]) ->
  # (metric > model > iteration > list(pred = , obs = ))
  transposed <- lapply(cfscore$metrics, function(m) {
    P <- lapply(names(cfscore$predictions), function(p) {
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
    names(P) <- names(cfscore$predictions)
    P
  })
  names(transposed) <- cfscore$metrics

  # # summarize
  conf.int <- lapply(cfscore$metrics, function(m) {
    CI <- lapply(names(cfscore$predictions), function(p) {
      if (m != "calplot") {
        return(ci(transposed[[m]][[p]], cover = 0.95))
      } else {
        return(NA)
      }
    })
    names(CI) <- names(cfscore$predictions)
    CI
  })
  names(conf.int) <- cfscore$metrics

  list(
    results = conf.int,
    raw = transposed
  )
}


lapply_progress <- function(x, FUN, task_description, progress = TRUE) {

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
