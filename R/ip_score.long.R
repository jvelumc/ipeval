#' @import data.table
ip_score.long <- function(predictions, data_long, outcome, treatment_formula,
                         treatment_of_interest,
                         metrics = c("auc", "brier", "oeratio", "calplot"),
                         time_horizon, cens.model = "cox",
                         null.model = TRUE, stable_iptw = FALSE,
                         bootstrap = 0, bootstrap_progress = TRUE,
                         iptw, ipcw, quiet = FALSE) {

  data_long <- copy(data_long)

  ipt <- ipt_weights(data_long, treatment_formula)

  data_long[, ipt.visit := ipt$weights]
  data_long[, iptw := cumprod(ipt.visit), by = id]
  data_long[, trt := all(A == treatment_of_interest[seq_len(.N)]), by = id]
  data_flat <- data_long[, .SD[.N], by = id]
  cfscore <- ip_score(predictions, data_flat, substitute(outcome),
                     treatment_formula = trt ~ 1, treatment_of_interest = 1,
                     metrics = metrics, time_horizon = time_horizon,
                     iptw = data_flat$iptw, ipcw = rep(1, nrow(data_flat)),
                     quiet = TRUE)
  return(cfscore)
}
