#' Performance in observed dataset This function exists only to demonstrate the
#' difference between 'normal' performance and counterfactual performance. It is
#' not user friendly and should not be relied on. It does not support cox models
#' out of the box.
#'
#' @param object One of the following three options to be validated:
#' \itemize{
#'   \item a numeric vector, corresponding to risk predictions
#'   \item a glm model
#'   \item a (named) list, with one or more of the previous 2 options, for
#'   validating and comparing multiple models at once.
#' }
#' @param data A data.frame containing the observed outcome.
#' @param outcome The outcome, to be evaluated within data. This could either be
#'   the name of a numeric column in data, or a Surv object for time-to-event
#'   data, e.g. Surv(time, status), if time and status are columns in data.
#' @param metrics A character vector specifying which performance metrics to be
#'   computed. Options are c("auc", "brier", "oeratio", "calplot").
#'
#' @returns Performance metrics in the observed dataset.
#' @export
observed_score <- function(object, data, outcome,
                    metrics = c("auc", "brier", "oeratio", "calplot")) {

  cfscore <- list()

  # get the observed outcome
  cfscore$outcome <- extract_outcome(data, substitute(outcome))
  if (inherits(cfscore$outcome, "Surv")) {
    stop("This function does not support survival data")
  } else {
    cfscore$outcome_type <- "binary"
  }

  # make a list of risk predictions
  object <- make_list_if_not_list(object)
  cfscore$predictions <- lapply(
    X = object,
    FUN = function(x) {
      if (is.numeric(x) && is.null(dim(x))) {
        x # user supplied risk predictions
      } else {
        stats::predict(x, newdata = data, type = "response")
      }
    }
  )

  cfscore$ipt$weights <- rep(1, nrow(data))
  cfscore$observed_treatment <- rep(1, nrow(data))
  cfscore$treatment_of_interest <- 1
  cfscore$metrics <- metrics
  cfscore$score <- get_metrics(cfscore)
  cfscore$bootstrap_iterations <- 0
  cfscore$quiet <- TRUE
  # rearrange such that score is the first entry of the list
  front <- c("score")
  cfscore <- cfscore[c(front, setdiff(names(cfscore), front))]
  class(cfscore) <- "CFscore"

  return(cfscore)

}
