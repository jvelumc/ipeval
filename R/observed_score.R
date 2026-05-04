#' Performance in observed dataset
#'
#' This function computes the performance of the predictions in the given data,
#' which may contain a mix of treated and untreated subjects. It exists only to
#' demonstrate the difference between 'normal' performance and counterfactual
#' performance. It is not user friendly and should not be relied on. It does not
#' support time-to-event data.
#'
#' @param object One of the following three options to be validated:
#' \itemize{
#'   \item a numeric vector, corresponding to risk predictions
#'   \item a glm model
#'   \item a (named) list, with one or more of the previous 2 options, for
#'   validating and comparing multiple models at once.
#' }
#' @param data A data.frame containing the observed outcome.
#' @param outcome The outcome, to be evaluated within data. This should be
#'   the name of a numeric column in data.
#' @param metrics A character vector specifying which performance metrics to be
#'   computed. Options are c("auc", "brier", "oeratio", "calplot").
#'
#' @returns Performance metrics in the observed dataset.
#' @export
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
#'
#' observed_score(
#'   object = list("ran" = random, "mod" = model),
#'   data = data,
#'   outcome = Y,
#'   metrics = c("auc", "brier", "oeratio")
#' )

observed_score <- function(object, data, outcome,
                    metrics = c("auc", "brier", "oeratio", "calplot")) {

  score <- list()

  # get the observed outcome
  score$outcome <- extract_outcome(data, substitute(outcome))
  if (inherits(score$outcome, "Surv")) {
    stop("This function does not support survival data")
  } else {
    score$outcome_type <- "binary"
  }

  # make a list of risk predictions
  object <- make_list_if_not_list(object)
  score$predictions <- lapply(
    X = object,
    FUN = function(x) {
      if (is.numeric(x) && is.null(dim(x))) {
        x # user supplied risk predictions
      } else {
        stats::predict(x, newdata = data, type = "response")
      }
    }
  )

  score$ipt$weights <- rep(1, nrow(data))
  score$observed_treatment <- rep(1, nrow(data))
  score$treatment_of_interest <- 1
  score$metrics <- metrics
  score$score <- get_metrics(score)
  score$bootstrap_iterations <- 0
  score$quiet <- TRUE
  # rearrange such that score is the first entry of the list
  front <- c("score")
  score <- score[c(front, setdiff(names(score), front))]
  class(score) <- "ip_score"

  return(score)

}
