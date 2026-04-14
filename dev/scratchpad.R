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


# Surv(time, status == 0) ~ 1
