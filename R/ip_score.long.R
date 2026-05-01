# assuming:
# - temperal ordering: P, L0, A0, Y1, L1, A1, Y2, ...
# - visit 0 is on time 0
# - visit times are same for everyone
# - visits are not necessarily equi-spaced?
# - treatment before time 0 is 0.
# - data_long has column 'id' indicating subject id's
# - data_long had column 'visit_id' indicating visit id, starting from 0

ip_score.long <- function(data_long, treatment_formula, treatment_of_interest,
                          n_visits) {



  ipt.visits <- ipt_weights(data_long, treatment_formula)
  ipt.product <- tapply(ipt.visits$weights, data_long$id,
                        FUN = prod)

  data_flat <- data.frame(
    id = unique(data_long$id),
    ipt = ipt.product,
    trt = faithful_to_trt(data_long, treatment_of_interest)
  )




  return(data_flat)
}

faithful_to_trt <- function(data_long, treatment_of_interest) {
  tapply(
    data_long$A == treatment_of_interest[data_long[["visit"]] + 1],
    data_long$id,
    FUN = all
  )
}



# ip_score.long <- function(predictions, data_long, outcome,
#                           treatment_formula, treatment_of_interest,
#                           metrics = c("auc", "brier", "oeratio", "calplot"),
#                           time_horizon, null.model = TRUE, id = id) {
#
#
#
# #   data_long <- copy(data_long)
# #
# #   ipt <- ipt_weights(data_long, treatment_formula)
# #
# #   data_long[, ipt.visit := ipt$weights]
# #   data_long[, iptw := cumprod(ipt.visit), by = id]
# #   data_long[, trt := all(A == treatment_of_interest[seq_len(.N)]), by = id]
# #   data_flat <- data_long[, .SD[.N], by = id]
# #   cfscore <- ip_score(predictions, data_flat, substitute(outcome),
# #                      treatment_formula = trt ~ 1, treatment_of_interest = 1,
# #                      metrics = metrics, time_horizon = time_horizon,
# #                      iptw = data_flat$iptw, ipcw = rep(1, nrow(data_flat)),
# #                      quiet = TRUE)
# #
#
#
#
#   ipt.visits <- ipt_weights(data_long, treatment_formula)
#   ipt.cumprod <- tapply(ipt.visit, )
#
#   construct_ip_object(
#     outcome = ,
#     treatment =,
#     predictions = ,
#     ipt =,
#     ipc = ,
#     metrics =
#   )
#
#   return(cfscore)
# }
