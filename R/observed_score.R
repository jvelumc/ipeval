#' Performance in observed dataset This function exists only to demonstrate the
#' difference between 'normal' performance and counterfactual performance. It is
#' not user friendly and should not be used.
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
        predict(x, newdata = data, type = "response")
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
