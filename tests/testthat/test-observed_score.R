test_that("observed_score works", {
  n <- 10000
  data <- data.frame(id = 1:n)
  data$L <- rnorm(n)
  data$A <- rbinom(n, 1, plogis(2*data$L))
  data$P <- rnorm(n)
  data$Z <- runif(n)
  data$Y0 <- rbinom(n, 1, plogis(0.5 + data$L + 1.25 * data$P))
  data$Y1 <- rbinom(n, 1, plogis(0.5 + data$L + 1.25 * data$P - 0.9*data$A))
  data$Y <- ifelse(data$A == 1, data$Y1, data$Y0)

  model <- glm(Y ~ L + A + P, family = "binomial", data = data)

  ipscore <- observed_score(model, data, Y, c("auc", "brier", "oeratio"))
  rrscore <- riskRegression::Score(list(model), Y ~ 1, data, null.model = FALSE)
  oeratio <- mean(data$Y)/mean(predict(model, type = "response"))

  expect_equal(ipscore$score$auc[[1]], rrscore$AUC$score$AUC)
  expect_equal(ipscore$score$brier[[1]], rrscore$Brier$score$Brier)
  expect_equal(ipscore$score$oeratio[[1]], oeratio)
})
