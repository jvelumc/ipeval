#' Counterfactual validation score
#'
#' Estimates the predictive performance of predictions under interventions, by
#' forming a weighted pseudopopulation in which every subject was assigned the
#' treatment of interest.
#'
#' When supplying a glm or coxph model, the model should be able to estimate
#' risks under the intervention of interest. This could be done in two ways: the
#' model does not have a treatment covariate, and always estimates the risk
#' under intervention of interest. Alternatively, the model has a covariate for
#' treatment. This function then automatically estimates the risk under the
#' treatment of interest for all subjects, even if they were assigned
#' alternative treatment.
#'
#' All performance metrics are computed on the weighted population in which
#' every subject was counterfactually assigned treatment of interest. \code{auc}
#' is area under the (ROC) curve. Brier score is defined as \code{1 / sum(iptw)
#' sum(predictions_i - outcome_i)^2}. Scaled brier score is also available
#' (\code{metrics = "scaled_brier"}). For the O/E ratio, the numerator
#' (observed) is the weighted fraction of events in the pseudopopulation, and
#' the denominator (expected) is the unweighted mean of risk estimates of the
#' original unweighted population. The \code{calplot} option generates a
#' calibration plot, with default 8 knots. More/less knots can be specified by
#' appending calplot with a number indicating the number of knots, e.g.
#' \code{metrics = calplot10} for 10 knots.
#'
#' Stabilized IPT-weigths are computed by estimating a null model for treatment.
#' E.g. weights are \code{P(A = a) / P(A = a | L = l)}, if the given
#' treatment_formula is \code{A ~ L}.
#'
#' @param object One of the following three options to be validated:
#' \itemize{
#'   \item a numeric vector, corresponding to risk predictions under
#'   intervention of interest.
#'   \item a glm or coxph model, capable of estimating risks under intervention
#'   of interest. See details.
#'   \item a (named) list, with one or more of the previous 2 options, for
#'   validating and comparing multiple models at once.
#' }
#' @param data A data.frame containing the observed outcome, assigned treatment,
#'   and necessary confounders for the validation of object.
#' @param outcome The outcome, to be evaluated within data. This could either be
#'   the name of a numeric/logical column in data, or a Surv object for
#'   time-to-event data, e.g. Surv(time, status), if time and status are columns
#'   in data.
#' @param treatment_formula A formula which identifies the treatment (left hand
#'   side) and the confounders (right hand side) in the data. E.g. A ~ L. The
#'   confounders are used to estimate the inverse probability of treatment
#'   weights (IPTW) model. The IPTW can also be specified themselves using the
#'   iptw argument, in which case the right hand side of this formula is
#'   ignored.
#' @param treatment_of_interest A treatment level for which the counterfactual
#'   perormance measures should be evaluated.
#' @param metrics A character vector specifying which performance metrics to be
#'   computed. Options are c("auc", "brier", "oeratio", "calplot"). See details.
#' @param time_horizon For time to event data, the prediction horizon of
#'   interest.
#' @param cens_model Model for estimating inverse probability of censored
#'   weights (IPCW). Methods currently implemented are Kaplan-Meier ("KM") or
#'   Cox ("cox"), both applied to the censored times. KM is only supported when
#'   the right hand side of cens_formula is 1.
#' @param cens_formula Formula for which the r.h.s. determines the censoring
#'   probabilities. I.e. ~ x1 + x2.
#' @param null_model If TRUE fit a risk prediction model which ignores the
#'   covariates and predicts the same value for all subjects. The model is
#'   fitted using the data in which all subjects are counterfactually assigned
#'   the treatment of interest (using the IPTW, as estimated using the
#'   treatment_formula or as given by the iptw argument). For time-to-event
#'   outcomes, the subjects are also counterfactually uncensored (using the
#'   IPCW, as estimated using the cens_formula, or as given by the ipcw
#'   argument).
#' @param stable_iptw if TRUE, estimate stabilized IPT-weights. See details.
#' @param bootstrap If this is an integer greater than 0, this indicates the
#'   number of bootstrap iterations, to compute 95\% confidence intervals around
#'   the performance metrics.
#' @param bootstrap_progress if set to TRUE, print a progress bar indicating the
#'   progress of bootstrap procedure.
#' @param iptw A numeric vector, containing the inverse probability of treatment
#'   weights. These are normally computed using the treatment_formula, but they
#'   can be specified directly via this argument. If specified via this
#'   argument, bootstrap is not possible.
#' @param ipcw A numeric vector, containing the inverse probability of censor
#'   weights. These are normally computed using the cens_formula, but they can
#'   be specified directly via this argument. If specified via this argument,
#'   bootstrap is not possible.
#' @param quiet If set to TRUE, don't print assumptions.
#'
#' @returns A list with performance metrics, of class 'ip_score", for which the
#' print and plot methods are implemented.
#'
#' @export
#'
#' @references Keogh RH, Van Geloven N. Prediction Under Interventions:
#'   Evaluation of Counterfactual Performance Using Longitudinal Observational
#'   Data. Epidemiology. 2024;35(3):329-339.
#'
#'   Boyer CB, Dahabreh IJ, Steingrimsson JA. Estimating and Evaluating
#'   Counterfactual Prediction Models. Statistics in Medicine.
#'   2025;44(23-24):e70287.
#'
#'   Pajouheshnia R, Peelen LM, Moons KGM, Reitsma JB, Groenwold RHH. Accounting
#'   for treatment use when validating a prognostic model: a simulation study.
#'   BMC Medical Research Methodology. 2017;17(1):103.
#'
#' @examples
#' n <- 1000
#'
#' data <- data.frame(L = rnorm(n), P = rnorm(n))
#' data$A <- rbinom(n, 1, plogis(data$L))
#' data$Y <- rbinom(n, 1, plogis(0.1 + 0.5*data$L + 0.7*data$P - 2*data$A))
#'
#' random <- runif(n, 0, 1)
#' model <- glm(Y ~ A + P, data = data, family = "binomial")
#' naive_perfect <- data$Y
#'
#' score <- ip_score(
#'   object = list("ran" = random, "mod" = model, "per" = naive_perfect),
#'   data = data,
#'   outcome = Y,
#'   treatment_formula = A ~ L,
#'   treatment_of_interest = 0,
#' )
#' print(score)
#' plot(score)

ip_score <- function(object, data, outcome, treatment_formula,
                    treatment_of_interest,
                    metrics = c("auc", "brier", "oeratio", "calplot"),
                    time_horizon, cens_model = "KM", cens_formula = ~ 1,
                    null_model = TRUE, stable_iptw = FALSE,
                    bootstrap = 0, bootstrap_progress = TRUE,
                    iptw, ipcw, quiet = FALSE) {


  # checking inputs ---------------------------------------------------------

  check_missing(object)
  check_missing(data)
  check_missing(outcome)
  check_missing(treatment_formula)
  check_missing(treatment_of_interest)

  if (cens_model == "KM")
    stopifnot("censoring model must be ~ 1 if modeled via KM" =
                rhs_is_one(cens_formula))


  cens_model <- match.arg(cens_model, choices = c("cox", "KM"))

  # assert rhs(outcome_formula != 1) iff surv model AND!missing(iptw_weights)
  # assert longest surv time is longer than time horizon, to avoid annoying
  # weights

  if (bootstrap != 0)
    stopifnot("can't bootstrap if iptw are given" = missing(iptw))


  score <- list()

  # get the observed outcome
  score$outcome <- extract_outcome(data, substitute(outcome))
  if (inherits(score$outcome, "Surv")) {
    score$outcome_type <- "survival"
    score$time_horizon <- time_horizon
    score$status_at_horizon <- ifelse(
      test = score$outcome[, 1] <= time_horizon,
      yes = score$outcome[, 2],
      no = FALSE
    )
  } else {
    score$outcome_type <- "binary"
  }

  # get the treatment
  score$treatment_column <- treatment_formula[[2]]
  score$observed_treatment <- extract_lhs(data, treatment_formula)
  score$treatment_of_interest <- treatment_of_interest
  stopifnot("Treatment is not binary" =
              setequal(unique(score$observed_treatment), c(0,1)))
  stopifnot("Treatment_of_interest must be either 0 or 1" =
    treatment_of_interest == 0 || treatment_of_interest == 1)

  # make a list of risk predictions
  object <- make_list_if_not_list(object)
  score$predictions <- lapply(
    X = object,
    FUN = function(x) {
      if (is.numeric(x) && is.null(dim(x))) {
        stopifnot("Predictions must be of length nrow(data)" =
                    length(x) == nrow(data))
        stopifnot("Predictions must be in interval [0,1]" =
                    all(x >= 0) && all(x <= 1))
        x # user supplied risk predictions
      } else {
        predict_CF(
          x,
          data,
          score$treatment_column,
          score$treatment_of_interest,
          time_horizon = score$time_horizon
        )
      }
    }
  )

  # get iptw
  score$ipt$method = "weights manually specified"
  if (missing(iptw)) {
    score$ipt$method <- "binomial glm"
    score$ipt$confounders <- all.vars(treatment_formula)[-1]
    score$ipt$propensity_formula <- treatment_formula
    ipt <- ipt_weights(data, treatment_formula)
    iptw <- ipt$weights
    score$ipt$model <- ipt$model

    if (stable_iptw == TRUE) {
      stable_treatment_formula <-
        stats::update.formula(treatment_formula, . ~ 1)
      sipt <- ipt_weights(data, stable_treatment_formula)
      iptw <- 1/sipt$weights * iptw
      score$ipt$stable_model <- sipt$model
    }
  }
  score$ipt$weights <- iptw

  # get ipcw
  if (score$outcome_type == "survival") {
    score$ipc$method <- "weights manually specified"
    if (missing(ipcw)) {
      score$ipc$method <- cens_model

      # annoying code to combine the Surv object from the outcome with the given
      # r.h.s. of the cens_formula
      score$ipc$cens_formula <- stats::update.formula(
        old = cens_formula,
        new = substitute(outcome ~ ., list(outcome = substitute(outcome)))
      )
      ipc <- ipc_weights(data, score$ipc$cens_formula,
                         cens_model, time_horizon)
      ipcw <- ipc$weights
      score$ipc$model <- ipc$model

    }
    score$ipc$weights <- ipcw
  }

  # add null model if required
  if (null_model == TRUE) {
    pseudo_ids <- score$observed_treatment == score$treatment_of_interest
    if (score$outcome_type == "binary") {
      null_model <- stats::lm(
        score$outcome[pseudo_ids] ~ 1,
        weights = score$ipt$weights[pseudo_ids]
      )
      null_preds <- stats::predict.lm(null_model, newdata = data)
    } else { # survival outcome
      # fit null model in a pseudopopuluation where everyone was treated and
      # no censoring
      # this seems wrong, but we only care about horizon, for which this works
      # (i think)
      uncensor_ids <- score$ipc$weights != 0
      cf_ids <- pseudo_ids & uncensor_ids

      null_model <- stats::weighted.mean(
        score$status_at_horizon[cf_ids],
        score$ipt$weights[cf_ids]*score$ipc$weights[cf_ids]
      )

      null_preds <- rep(null_model, nrow(data))
    }
    score$predictions <- c(
      list("null model" = null_preds),
      score$predictions
    )
  }


  score$metrics <- metrics
  score$score <- get_metrics(score)

  score$bootstrap_iterations <- bootstrap
  if (bootstrap != 0) {
    if (score$ipt$method == "weights manually specified" ||
        (score$outcome_type == "survival" &&
         score$ipc$method == "weights manually specified")) {
      stop("Can not compute confidence intervals if manually specifying
               iptw/ipcw.")
    }
    stopifnot(
      "bootstrap should be the number of iterations and must be greater than 1" =
        bootstrap > 1)

    score$bootstrap_progress <- bootstrap_progress
    score$bootstrap <- bootstrap(data, score)
  }

  score$quiet <- quiet


  # rearrange such that score is the first entry of the list
  front <- c("score")
  score <- score[c(front, setdiff(names(score), front))]
  class(score) <- "ip_score"

  return(score)
}

get_metrics <- function(ip_score) {
  metrics <- list()

  if (ip_score$outcome_type == "survival") {
    weights <- ip_score$ipt$weights * ip_score$ipc$weights
    outcome <- ip_score$status_at_horizon
  } else {
    weights <- ip_score$ipt$weights
    outcome <- ip_score$outcome
  }

  for (m in ip_score$metrics) {
    metrics[[m]] <- sapply(
      X = ip_score$predictions,
      FUN = function(x) {
        cf_metric(
          m,
          obs_outcome = outcome,
          obs_trt = ip_score$observed_treatment,
          cf_pred = x,
          cf_trt = ip_score$treatment_of_interest,
          ipw = weights,
          nullpred = ip_score$predictions[["null model"]]
        )
      }
    )
  }

  metrics
}

name_unnamed_list <- function(x) {
  # give names, if not named
  sapply(
    1:length(x),
    function(i)

      if (is.null(names(x)[i]) || names(x)[i] == "") {
        paste0("model.", i)
      } else {
        names(x[i])
      }
  )
}

make_list_if_not_list <- function(x) {
  if (!("list" %in% class(x)))
    x <- list(x)
  names(x) <- name_unnamed_list(x)
  x
}

extract_outcome <- function(data, outcome) {
  # attempt to extract the outcome from the data, and perform various sanity
  # checks

  y <- tryCatch(
    eval(outcome, envir = data),
    error = function(e) {
      stop(sprintf("Outcome %s not found in data", deparse(outcome)),
           call. = FALSE)
    }
  )

  if (!( ((is.numeric(y) || is.logical(y)) && is.vector(y)) ||
         inherits(y, "Surv") )) {
    stop("Outcome must be a numeric vector or a Surv object", call. = FALSE)
  }

  if (length(y) != nrow(data)) {
    stop("Outcome must be of length nrow(data)", call. = FALSE)
  }

  if (inherits(y, "Surv")) {
    time  <- y[, 1]
    event <- y[, 2]

    if (any(!is.finite(time)) || any(time < 0)) {
      stop("Survival times must be finite & nonnegative", call. = FALSE)
    }

    if (!all(event %in% c(0, 1))) {
      stop("Event indicator must be binary (0/1)", call. = FALSE)
    }
  } else {
    if (!all(y %in% c(0, 1))) {
      stop("Outcome must be binary (0/1)", call. = FALSE)
    }
  }

  y
}
