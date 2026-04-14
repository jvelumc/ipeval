data <- data.frame(
  time = c(1,2,2,3),
  status = c(1, 1, 0, 1)
)

# force censors after events by increasing timing of censors by epsilon
data2 <- data.frame(
  time = c(1,2,2.0000001,3),
  status = c(1, 1, 0, 1)
)


get_step <- function(fit) {
  stepfun(fit$time, c(1, fit$surv), right = TRUE)
}

# 'truth'
get_step(survfit(Surv(time, !status) ~ 1, data2))


get_step(survfit(Surv(time, !status) ~ 1, data))
get_step(prodlim(Surv(time, status) ~ 1, data, reverse = TRUE))
get_step(survfit(coxph(Surv(time, !status) ~ 1, data)))
get_step(survfit(coxph(Surv(time, !status) ~ 1, data2)))
# why is this not 0.5?
