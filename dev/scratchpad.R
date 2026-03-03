n <- 10000

data <- data.frame(id = 1:n)
data$L <- rnorm(n)
data$A <- rbinom(n, 1, plogis(2*data$L))
data$P <- rnorm(n)
data$Y0 <- rbinom(n, 1, plogis(0.5 + data$L + 1.25 * data$P))
data$Y1 <- rbinom(n, 1, plogis(0.5 + data$L + 1.25 * data$P - 0.9*data$A))
data$Y <- ifelse(data$A == 1, data$Y1, data$Y0)

trt_model <- glm(A ~ L, family = "binomial", data = data)
propensity_score <- predict(trt_model, type = "response")
data$iptw <- 1 / ifelse(data$A == 1, propensity_score, 1 - propensity_score)
causal_model <- glm(Y ~ A + P, family = "binomial", data = data, weights = iptw)

CFscore(
  object = causal_model,
  data = data,
  outcome_formula = Y ~ 1,
  treatment_formula = A ~ L,
  treatment_of_interest = 0,
  metrics = c("oeratio", "brier", "brier")
)

CFscore(
  object = causal_model,
  data = data,
  outcome_formula = Y ~ 1,
  treatment_formula = A ~ L,
  treatment_of_interest = 1,
  metrics = c("oeratio"),
  null.model = FALSE
)

observed_score(
  predict_CF(causal_model, data, "A", 1),
  data, Y1 ~ 1, metrics = "oeratio")


data$A <- 1

riskRegression::Score(
  object = list(causal_model),
  formula = Y1 ~ 1,
  data = data
)
