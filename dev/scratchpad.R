library(prodlim)
library(riskRegression)

n <- 10000

data <- data.frame(id = 1:n)
data$L <- rnorm(n)
data$A <- rbinom(n, 1, plogis(2*data$L))
data$P <- rnorm(n)
data$time_event1 <- simulate_time_to_event(n, 0.1, 0.5*data$P)
data$time_censor2 <- simulate_time_to_event(n, 0.1, 0.3*data$P)
data$time_censor3 <- ifelse(rbinom(n, 1, 0.5), 5, Inf)

m_times <- as.matrix(data[, c("time_event1", "time_censor2", "time_censor3")])

data$status <- max.col(-m_times)
data$time <- pmin(data$time_event1, data$time_censor2, data$time_censor3)

prodlim(Hist(time, status, cens.code = 3) ~ 1, data, reverse = TRUE)




my_func <- function(data, outcome, cens_formula) {

  # proposed approach:
  # if user wants 1 censoring mechanism, we assume event = 1, censor = 0.
  # The user then specifies cens_formula = ~ x1 + x2.
  # this is combined with outcome on the l.h.s.

  # if user wants more than 1 censoring mechanism, user must specify code for
  # each event/censoring mechanism.
  # The user specifies cens_formula then as a list as follows
  # cens_formula = list(
  #   Hist(time, status, cens.code = 2) ~ x1,
  #   Hist(time, status, cens.code = 3) ~ 1
  # )
  # and a method for each censoring mechanism i.e. c("cox", "km")

  cens_formula <- make_list_if_not_list(cens_formula)

  if (length(cens_formula) == 1) {
    updated_formula <- stats::update.formula(
      old = cens_formula[[1]],
      new = substitute(outcome ~ ., list(outcome = substitute(outcome)))
    )
    return(ipc_weights(data, updated_formula, "KM", 6))
  } else {
    lapply(cens_formula, FUN = function(f, data) {
      # run some adapted version of ipc_weights for each formula
    })
    # multiply the weights for every subject
  }
}

data$ipcw_onemech <- my_func(data, Surv(time, status == 1), cens_formula = ~ 1)$weights
