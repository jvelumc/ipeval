library(data.table)
library(survival)
# simulate some longitudinal trt time to event data

n <- 3000
n_visits <- 5

gamma_0 <- -1
gamma_L <- 0.5

alpha_0 <- -2
alpha_A <- -0.5
alpha_L <- 0.5
alpha_U <- 0.5

simulate_longitudinal <- function(n, fix_trt = NULL) {
  U <- rnorm(n, 0, 0.1)

  A <- matrix(nrow = n, ncol = n_visits)
  L <- matrix(nrow = n, ncol = n_visits)

  time <- rep(NA, n)
  status <- rep(NA, n)

  simulate_A <- function(i, L, fix_trt) {
    if (is.null(fix_trt)) {
      return(rbinom(n, 1, plogis(gamma_0 + gamma_L * L[, i])))
    } else {
      return(rep(fix_trt, n))
    }
  }

  L[, 1] <- rnorm(n, U, 1)
  A[, 1] <- simulate_A(1, L, fix_trt)

  for (i in 2:n_visits) {
    L[, i] <- rnorm(n, 0.8 * L[, i - 1] - A[, i - 1] + 0.1 * (i-i) + U)
    A[, i] <- ifelse(A[, i - 1] == 1, rep(1, n), simulate_A(i, L, fix_trt))
  }

  for (i in 1:n_visits){
    new.t <- simulate_time_to_event(
      n = n,
      constant_baseline_haz = 1,
      LP = alpha_0 + alpha_A * A[, i] + alpha_L * L[, i] + alpha_U * U
    )
    time <- ifelse(is.na(time) & new.t < 1, i - 1 + new.t,time)
  }
  status <- ifelse(is.na(time), 0, 1)
  time <- ifelse(is.na(time), 5, time)

  colnames(A) <- paste0("A", 0:(n_visits - 1))
  colnames(L) <- paste0("L", 0:(n_visits - 1))

  data.table(id = 1:n, time, status, A, L, U)
}

make_long <- function(data) {
  long <- melt(data,
               measure.vars = patterns("^A", "^L"),
               variable.name = "visit",
               value.name = c("A", "L"),
               variable.factor = F)

  long[, visit := as.numeric(visit) - 1]

  long <- long[order(id, visit)]
  long[, `:=`(
    time_start = visit,
    time_end = pmin(time, visit + 1)
  )]
  long <- long[time_start < time, ]
  long[, status := fifelse(time_end == time, status, 0)]
  long[, time := NULL]

  long[, paste0("A_lag", 1:(n_visits-1)) := shift(A, n = 1:(n_visits-1), fill = 0), by = .(id)]
  long[, paste0("L_lag", 1:(n_visits-1)) := shift(L, n = 1:(n_visits-1), fill = 0), by = .(id)]
  long[, L0 := first(L), by = .(id)]

  return(long[])
}

fit_model <- function(data_long) {
  data_long <- copy(data_long)
  ipt <- ipt_weights(data_long, A ~ L * A_lag1)
  data_long[, wt := ipt$weights]
  data_long[, ipw := cumprod(wt), by = id]

  coxph(
    Surv(time_start, time_end, status) ~ A + A_lag1 + A_lag2 + A_lag3 + A_lag4 + L0,
    data = data_long,
    weights = ipw
  )
}


simulate_iteration <- function() {
  df_dev <- simulate_longitudinal(n)
  df_val <- simulate_longitudinal(n)

  df_cf0 <- simulate_longitudinal(n, 0)
  df_cf1 <- simulate_longitudinal(n, 1)

  df_dev_long <- make_long(df_dev)

  # model dev ---------------------------------------------------------------

  cox.msm <- fit_model(df_dev_long)

  cumhaz <- basehaz(cox.msm, centered = FALSE)$hazard
  event.times <- basehaz(cox.msm, centered = FALSE)$time
  cumhaz.fun <- stepfun(event.times, c(0, cumhaz))

  # hazard rate term
  hrt <- function(variable, value) {
    exp(coef(cox.msm)[variable]*value)
  }

  # # risk under never treated
  risk0 <- function(t, L0) {
    1 - exp(-cumhaz.fun(t)*hrt("L0", L0))
  }

  # # risk under always treated
  risk1 <- function(t, L0) {
    1 - exp(-(
      cumhaz.fun(min(t, 1))*
        hrt("L0", L0)*hrt("A", 1) +

        (t >= 1) *
        (cumhaz.fun(min(t, 2)) - cumhaz.fun(1))*
        hrt("L0", L0)*hrt("A", 1)*hrt("A_lag1", 1) +

        (t >= 2) *
        (cumhaz.fun(min(t, 3)) - cumhaz.fun(2))*
        hrt("L0", L0)*hrt("A", 1)*hrt("A_lag1", 1)*hrt("A_lag2", 1) +

        (t >= 3) *
        (cumhaz.fun(min(t, 4)) - cumhaz.fun(3))*
        hrt("L0", L0)*hrt("A", 1)*hrt("A_lag1", 1)*hrt("A_lag2", 1)*hrt("A_lag3", 1) +

        (t >= 4) *
        (cumhaz.fun(min(t, 5)) - cumhaz.fun(4))*
        hrt("L0", L0)*hrt("A", 1)*hrt("A_lag1", 1)*hrt("A_lag2", 1)*
        hrt("A_lag3", 1)*hrt("A_lag4", 1)
    ))
  }

  # validate
  risk_under_0 <- risk0(5, df_val$L0)
  risk_under_1 <- risk1(5, df_val$L0)
  df_val_long <- make_long(df_val)

  cfscore0 <- ip_score.long(predictions = risk_under_0, data_long = df_val_long,
                           outcome = Surv(time_end, status),
                           treatment_formula = A ~ A_lag1 * L,
                           treatment_of_interest = c(0,0,0,0,0),
                           metrics = c("oeratio", "auc", "brier"),
                           time_horizon = 5)

  true0 <- observed_score(risk0(5, df_cf0$L0), data = df_cf0, outcome = status,
                          metrics = c("oeratio", "auc", "brier"))


  cfscore1 <- ip_score.long(predictions = risk_under_1, data_long = df_val_long,
                           outcome = Surv(time_end, status),
                           treatment_formula = A ~ A_lag1 * L,
                           treatment_of_interest = c(1,1,1,1,1),
                           metrics = c("oeratio", "auc", "brier"),
                           time_horizon = 5)

  true1 <- observed_score(risk1(5, df_cf1$L0), data = df_cf1, outcome = status,
                          metrics = c("oeratio", "auc", "brier"))

  # this is crazy
  as.data.table(as.list(unlist(list(truth0 = sapply(true0$score, function(x) x[[1]]),
                                    cf0 = sapply(cfscore0$score, function(x) x[[2]]),
                                    truth1 = sapply(true1$score, function(x) x[[1]]),
                                    cf1 = sapply(cfscore1$score, function(x) x[[2]])
  ))))
}

results <- lapply_progress(1:200, function(x) simulate_iteration(), "simulating cox 1")
results <- rbindlist(results)

metrics <- c("oeratio", "auc", "brier")
for (m in metrics) {
  results[, paste0("bias0.", m) :=
            get(paste0("truth0.", m)) - get(paste0("cf0.", m))]

  results[, paste0("bias1.", m) :=
            get(paste0("truth1.", m)) - get(paste0("cf1.", m))]
}
results[, lapply(.SD, mean)]
