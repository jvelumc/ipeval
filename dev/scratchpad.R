n <- 10000

data <- data.frame(id = 1:n)
data$L <- rnorm(n)
data$A <- 1 + rbinom(n, 1, plogis(2*data$L))
data$P <- rnorm(n)
data$Z <- runif(n)
data$Y0 <- rbinom(n, 1, plogis(0.5 + data$L + 1.25 * data$P))
data$Y1 <- rbinom(n, 1, plogis(0.5 + data$L + 1.25 * data$P - 0.9*data$A))
data$Y <- ifelse(data$A == 2, data$Y1, data$Y0)

trt_model <- glm(A ~ L, family = "binomial", data = data)
propensity_score <- predict(trt_model, type = "response")
data$iptw <- 1 / ifelse(data$A == 1, propensity_score, 1 - propensity_score)
causal_model <- glm(Y ~ A + P, family = "binomial", data = data, weights = iptw)






ip_score(causal_model, data, Y, A ~ L, 0, metrics = "calplot8")

