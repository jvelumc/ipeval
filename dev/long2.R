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

df_dev <- simulate_longitudinal(n)


# input 1: outcome data
outcome_data <- df_dev[, c("id", "time", "status")]

# input 2: treatment/confounder data
# some wide to long function
df_wide <- wide_to_long(df_dev,
             baseline_variables = c("id", "U"),
             wide_variables = list(A = c("A0", "A1", "A2", "A3", "A4"),
                                   L = c("L0", "L1", "L2", "L3", "L4")),
             visit_times = c(0,1,2,3,4),
             outcome_times = outcome_data$time)

df_wide <- add_lag_terms(df_wide, "A", 1)
df_wide <- add_lag_terms(df_wide, "L", 1)

