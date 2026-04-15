# input checks

test_that("wrong input throws sensible errors", {
  n <- 1000
  adminstrative_censor <- 10
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)
  data$status <- ifelse(data$time_uncensored <= adminstrative_censor, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE, data$time_uncensored, adminstrative_censor)

  predictions <- runif(n, 0, 1)

  # object ------------------------------------------------------------------
  expect_error(
    CFscore(predictions[1:999], data = data, outcome = status, A ~ L, 1),
    "Predictions must be of length nrow"
  )

  expect_error(
    CFscore(lm(status ~ A, data), data = data, outcome = status, A ~ L, 1),
    "model class lm not supported"
  )

  expect_error(
    CFscore(runif(n, -1, 2), data = data, outcome = status, A ~ L, 1),
    "Predictions must be in interval"
  )

  # outcome -----------------------------------------------------------------
  expect_error(
    CFscore(predictions, data = data, outcome = B, A ~ L, 1),
    "Outcome B not found"
  )

  expect_error(
    CFscore(predictions, data = data, outcome = data$status[1:999], A ~ L, 1),
    "Outcome must be of length"
  )

  expect_error(
    CFscore(predictions, data = data, outcome = time, A ~ L, 1),
    "Outcome must be binary"
  )

  expect_error(
    CFscore(predictions, data = data, outcome = rep(1, 1000), A ~ L, 1,
            metrics = "auc"),
    "no controls"
  )

  # treatment
  expect_error(
    CFscore(predictions, data, status, A + B ~ L, 1),
    "treatment formula must be one variable"
  )

  expect_error(
    CFscore(predictions, data, status, ~ L, 1),
    "treatment formula must be one variable"
  )

  expect_error(
    CFscore(predictions, data, status, time ~ L, 1),
    "Treatment is not binary"
  )

  expect_error(
    CFscore(predictions, data, status, A ~ L, 2),
    "Treatment_of_interest"
  )
})

test_that("supplying (list of) model or predictions equivalent", {
  set.seed(1)
  n <- 1000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.5+0.2*data$L))
  data$Y <- rbinom(n, 1, plogis(0.3*data$L + 0.6*data$P - 0.5*data$A))

  model1 <- glm(Y ~ P, family = "binomial", data = data)
  model2 <- glm(Y ~ A + P, family = "binomial", data = data)

  expect_equal(
    CFscore(
      data = data,
      object = model1,
      outcome = Y,
      treatment_formula = A ~ L,
      treatment_of_interest = 0,
    ),
    CFscore(
      data = data,
      object = list(model1),
      outcome = Y,
      treatment_formula = A ~ L,
      treatment_of_interest = 0
    )
  )

  expect_equal(
    CFscore(
      data = data,
      object = list(model2),
      outcome = Y,
      treatment_formula = A ~ L,
      treatment_of_interest = 0
    ),
    CFscore(
      data = data,
      object = predict_CF(model2, data, "A", 0),
      outcome = Y,
      treatment_formula = A ~ L,
      treatment_of_interest = 0
    )
  )

  expect_equal(
    CFscore(
      data = data,
      object = list(aa = model2),
      outcome = Y,
      treatment_formula = A ~ L,
      treatment_of_interest = 0
    ),
    CFscore(
      data = data,
      object = list(aa = predict_CF(model2, data, "A", 0)),
      outcome = Y,
      treatment_formula = A ~ L,
      treatment_of_interest = 0
    )
  )
})


# TODO manual iptw/ipcw equivalent to equivalent models



# metrics

test_that("CFscore metrics equal to unobserved CF metrics, binary outcome", {

  set.seed(1)
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2+0.5*data$L))
  data$Y0 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P))
  data$Y1 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P - 0.4))
  data$Y <- ifelse(data$A == 1, data$Y1, data$Y0)

  model <- suppressWarnings(
    glm(
      Y ~ A + P,
      family = "binomial",
      data = data,
      weights = ipt_weights(data, A ~ L)$weights
    )
  )

  cfscore <- CFscore(
    data = data,
    object = model,
    outcome = Y,
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    null_model = FALSE
  )

  Y0_predicted <- predict_CF(model, data, "A", 0)
  score <- riskRegression::Score(
    list(Y0_predicted),
    formula = Y0 ~ 1,
    data = data,
    null.model = F
  )
  score$oe <- mean(data$Y0)/mean(Y0_predicted)

  expect_equal(unname(cfscore$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(cfscore$score$brier), score$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(cfscore$score$oeratio), score$oe, tolerance = 0.01)
})

test_that("CFscore metrics equal to unobserved CF metrics, surv, uncensored", {
  set.seed(1)
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))
  data$status <- 1 # no censoring, so status is always 1 at the end

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$time <- ifelse(data$A == 1, data$time1, data$time0)

  summary(data$time0)
  summary(data$time1)

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data,
    weights = ipt_weights(data, A ~ L)$weights
  )

  horizon <- 10

  cfscore <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = horizon,
    cens_model = "KM",
    cens_formula = ~ 1,
    null_model = FALSE
  )

  time0_predicted <- predict_CF(model, data, "A", 0, horizon)
  score <- riskRegression::Score(
    list(time0_predicted),
    formula = Hist(time0, status) ~ 1,
    data = data,
    null.model = F,
    times = horizon
  )
  score$oe <- mean(data$time0 <= horizon)/mean(time0_predicted)

  expect_equal(unname(cfscore$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(cfscore$score$brier), score$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(cfscore$score$oeratio), score$oe, tolerance = 0.01)
})

test_that("CFscore metrics equal to unobserved CF metrics, surv, censor at T", {
  set.seed(1)
  horizon <- 9.999
  adminstrative_censor <- 10
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)

  summary(data$time0)
  summary(data$time1)

  data$status <- ifelse(data$time_uncensored <= adminstrative_censor, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE, data$time_uncensored, adminstrative_censor)

  data$status_uncensored <- 1

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data
  ) # naive model, i.e. not adjusting for confounding

  cfscore_km <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = horizon,
    cens_model = "KM",
    cens_formula = ~ 1,
    null_model = FALSE
  )

  cfscore_cox <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = horizon,
    cens_model = "cox",
    cens_formula = ~ 1,
    null_model = FALSE
  )

  expect_equal(
    cfscore_km$ipc$weights,
    rep(1, n)
  )
  expect_equal(
    cfscore_cox$ipc$weights,
    rep(1, n)
  )

  time0_predicted <- predict_CF(model, data, "A", 0, horizon)
  score <- riskRegression::Score(
    list(time0_predicted),
    formula = Hist(time0, status_uncensored) ~ 1,
    data = data,
    null.model = F,
    times = horizon
  )
  score$oe <- mean(data$time0 <= horizon)/mean(time0_predicted)

  expect_equal(unname(cfscore_km$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(cfscore_km$score$brier), score$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(cfscore_km$score$oeratio), score$oe, tolerance = 0.01)

  expect_equal(unname(cfscore_cox$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(cfscore_cox$score$brier), score$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(cfscore_cox$score$oeratio), score$oe, tolerance = 0.01)
})

test_that("CFscore metrics equal to unobserved CF metrics, surv, KM censor", {
  set.seed(1)
  horizon <- 10
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$censortime <- simulate_time_to_event(n, 0.04, 0)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)
  data$status_uncensored <- 1

  summary(data$time0)
  summary(data$time1)
  summary(data$censortime)

  data$status <- ifelse(data$time_uncensored <= data$censortime, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE,
                      data$time_uncensored,
                      data$censortime)

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data
  )

  cfscore <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = horizon,
    cens_model = "KM",
    cens_formula = ~ 1,
    null_model = FALSE
  )

  time0_predicted <- predict_CF(model, data, "A", 0, horizon)
  score <- riskRegression::Score(
    list(time0_predicted),
    formula = Hist(time0, status_uncensored) ~ 1,
    data = data,
    null.model = F,
    times = horizon
  )
  score$oe <- mean(data$time0 <= horizon)/mean(time0_predicted)

  expect_equal(unname(cfscore$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(cfscore$score$brier), score$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(cfscore$score$oeratio), score$oe, tolerance = 0.01)

  # also try treatment == 1 for good measure
  cfscore1 <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 1,
    time_horizon = horizon,
    cens_model = "KM",
    cens_formula = ~ 1,
    null_model = FALSE
  )

  time1_predicted <- predict_CF(model, data, "A", 1, horizon)
  score1 <- riskRegression::Score(
    list(time1_predicted),
    formula = Hist(time1, status_uncensored) ~ 1,
    data = data,
    null.model = F,
    times = horizon
  )
  score1$oe <- mean(data$time1 <= horizon)/mean(time1_predicted)

  expect_equal(unname(cfscore1$score$auc), score1$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(cfscore1$score$brier), score1$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(cfscore1$score$oeratio), score1$oe, tolerance = 0.01)
})

test_that("CFscore metrics equal to unobserved CF metrics, surv, cox censor", {
  set.seed(1)
  horizon <- 10
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$censortime <- simulate_time_to_event(n, 0.04, 0.5*data$L + 0.6*data$P +
                                              0.2*data$A)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)
  data$status_uncensored <- 1

  summary(data$time0)
  summary(data$time1)
  summary(data$censortime)

  data$status <- ifelse(data$time_uncensored <= data$censortime, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE,
                      data$time_uncensored,
                      data$censortime)

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data
  )

  cfscore <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = horizon,
    cens_model = "cox",
    cens_formula = ~ L + P + A,
    null_model = FALSE
  )

  time0_predicted <- predict_CF(model, data, "A", 0, horizon)
  score <- riskRegression::Score(
    list(time0_predicted),
    formula = Hist(time0, status_uncensored) ~ 1,
    data = data,
    null.model = F,
    times = horizon
  )
  score$oe <- mean(data$time0 <= horizon)/mean(time0_predicted)

  expect_equal(unname(cfscore$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(cfscore$score$brier), score$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(cfscore$score$oeratio), score$oe, tolerance = 0.01)
})

test_that("CFscore metrics equal to unobserved CF metrics, surv, KM censor, stable weights", {
  set.seed(1)
  horizon <- 10
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$censortime <- simulate_time_to_event(n, 0.04, 0)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)
  data$status_uncensored <- 1

  data$status <- ifelse(data$time_uncensored <= data$censortime, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE,
                      data$time_uncensored,
                      data$censortime)

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data
  )

  cfscore <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = horizon,
    cens_model = "KM",
    cens_formula = ~ 1,
    stable_iptw = TRUE,
    null_model = FALSE
  )

  stable_min <- min(cfscore$ipt$weights)
  stable_max <- max(cfscore$ipt$weights)

  cfscore_unstable <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = horizon,
    cens_model = "KM",
    cens_formula = ~ 1,
    stable_iptw = FALSE,
    null_model = FALSE
  )

  unstable_min <- min(cfscore_unstable$ipt$weights)
  unstable_max <- max(cfscore_unstable$ipt$weights)

  expect_true(stable_min <= unstable_min)
  expect_true(stable_max <= unstable_max)


  time0_predicted <- predict_CF(model, data, "A", 0, horizon)
  score <- riskRegression::Score(
    list(time0_predicted),
    formula = Hist(time0, status_uncensored) ~ 1,
    data = data,
    null.model = F,
    times = horizon
  )
  score$oe <- mean(data$time0 <= horizon)/mean(time0_predicted)

  expect_equal(unname(cfscore$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(cfscore$score$brier), score$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(cfscore$score$oeratio), score$oe, tolerance = 0.01)

})


# minor bootstrap tests
test_that("results are in between lower & upper bootstrap", {
  set.seed(1)
  n <- 1000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2+0.5*data$L))
  data$Y0 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P))
  data$Y1 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P - 0.4))
  data$Y <- ifelse(data$A == 1, data$Y1, data$Y0)

  model <- suppressWarnings(
    glm(
      Y ~ A + P,
      family = "binomial",
      data = data
    )
  )

  cfscore <- CFscore(
    data = data,
    object = model,
    outcome = Y,
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    bootstrap = 200,
    null_model = FALSE
  )

  expect_true(
    cfscore$score$auc > cfscore$bootstrap$results$auc[[1]][1] &
      cfscore$score$auc < cfscore$bootstrap$results$auc[[1]][2]
  )
  expect_true(
    cfscore$score$brier > cfscore$bootstrap$results$brier[[1]][1] &
      cfscore$score$brier < cfscore$bootstrap$results$brier[[1]][2]
  )
  expect_true(
    cfscore$score$oeratio > cfscore$bootstrap$results$oeratio[[1]][1] &
      cfscore$score$oeratio < cfscore$bootstrap$results$oeratio[[1]][2]
  )
})

test_that("results are in between lower & upper bootstrap, surv, cox censor", {
  set.seed(1)
  horizon <- 10
  n <- 10000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$censortime <- simulate_time_to_event(n, 0.04, 0.5*data$L + 0.6*data$P +
                                              0.2*data$A)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)
  data$status_uncensored <- 1

  data$status <- ifelse(data$time_uncensored <= data$censortime, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE,
                      data$time_uncensored,
                      data$censortime)

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data
  )

  cfscore <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    cens_model = "cox",
    time_horizon = horizon,
    cens_formula = ~ L + P + A,
    bootstrap = 100,
    null_model = FALSE
  )

  expect_true(
    cfscore$score$auc > cfscore$bootstrap$results$auc[[1]][1] &
      cfscore$score$auc < cfscore$bootstrap$results$auc[[1]][2]
  )
  expect_true(
    cfscore$score$brier > cfscore$bootstrap$results$brier[[1]][1] &
      cfscore$score$brier < cfscore$bootstrap$results$brier[[1]][2]
  )
  expect_true(
    cfscore$score$oeratio > cfscore$bootstrap$results$oeratio[[1]][1] &
      cfscore$score$oeratio < cfscore$bootstrap$results$oeratio[[1]][2]
  )
})

# null model

test_that("null model binary outcome", {
  set.seed(1)
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2+0.5*data$L))
  data$Y0 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P))
  data$Y1 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P - 0.4))
  data$Y <- ifelse(data$A == 1, data$Y1, data$Y0)

  model <- suppressWarnings(
    glm(
      Y ~ A + P,
      family = "binomial",
      data = data
    )
  )

  nullmodel <- glm(Y0 ~ 1, data = data)

  cfscore <- CFscore(model, data, Y, A ~ L, 0,
                     metrics = c("brier", "scaled_brier"), null_model = TRUE)
  expect_equal(
    unname(nullmodel$coefficients[1]),
    cfscore$predictions$`null model`[[1]],
    tolerance = 0.01
  )

  # scaled brier score of null model is 0
  expect_equal(
    cfscore$score$scaled_brier[[1]],
    0
  )

  # scaled brier score of model is correct
  expect_equal(
    (1 - cfscore$score$brier[[2]]/cfscore$score$brier[[1]])*100,
    cfscore$score$scaled_brier[[2]]
  )
})

test_that("null model survival outcome", {
  set.seed(2)
  horizon <- 10
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$censortime <- simulate_time_to_event(n, 0.04, 0.5*data$L + 0.6*data$P +
                                              0.2*data$A)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)
  data$status_uncensored <- 1

  data$status <- ifelse(data$time_uncensored <= data$censortime, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE,
                      data$time_uncensored,
                      data$censortime)


  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    x = TRUE,
    data = data
  )
  nullmodel <- survival::survfit(survival::Surv(time0, status_uncensored) ~ 1, data)

  cfscore <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    metrics = c(),
    cens_model = "cox",
    cens_formula = ~ L + P + A,
    time_horizon = horizon,
    null_model = TRUE
  )

  expect_equal(
    cfscore$predictions$`null model`[[1]],
    1 - summary(nullmodel, times = horizon)$surv,
    tolerance = 0.01
  )

})

# ties

test_that("CFscore handles horizon on censor time correctly", {
  set.seed(1)

  adminstrative_censor <- 10
  n <- 10000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)

  data$status <- ifelse(data$time_uncensored <= adminstrative_censor, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE, data$time_uncensored, adminstrative_censor)

  data$status_uncensored <- 1

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data
  )

  cfscore_km999 <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = 9.999,
    cens_model = "KM",
    cens_formula = ~ 1,
    null_model = TRUE
  )

  cfscore_km10 <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = 10,
    cens_model = "KM",
    cens_formula = ~ 1,
    null_model = TRUE
  )

  cfscore_cox999 <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = 9.999,
    cens_model = "cox",
    cens_formula = ~ 1,
    null_model = TRUE
  )

  cfscore_cox10 <- CFscore(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = 10,
    cens_model = "cox",
    cens_formula = ~ 1,
    null_model = TRUE
  )

  expect_equal(cfscore_km999$score, cfscore_km10$score)
  expect_equal(cfscore_cox999$score, cfscore_cox10$score)
})





