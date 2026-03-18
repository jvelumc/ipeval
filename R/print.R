# pretty print, respecting output width, adding \n at the end
pp <- function(...) {
  txt <- paste0(c(...), collapse = "")
  cat(paste0(strwrap(txt), "\n"))
}


#' @export
print.CFscore <- function(x, ...) {
  if (x$quiet != TRUE) {
    assumptions(x)
  }
  numeric_metrics <- x$metrics[x$metrics != "calplot"]

  if (x$bootstrap_iterations != 0) {
    for (metric in numeric_metrics) {
      cat("\n", metric, "\n\n", sep = "")
      tab <- data.frame(model = names(x$predictions))
      tab[[metric]] <- x$score[[metric]]
      tab$lower <- sapply(x$bootstrap$results[[metric]], function(x) x[[1]])
      tab$upper <- sapply(x$bootstrap$results[[metric]], function(x) x[[2]])
      print(tab, digits = 3, row.names = FALSE)
    }
  } else {
    cat("\n")
    tab <- data.frame(model = names(x$predictions))
    for (metric in numeric_metrics) {
      tab[[metric]] <- x$score[[metric]]
    }
    print(tab, digits = 3, row.names = FALSE)
  }

  if ("calplot" %in% x$metrics) {
    plot(x, ...)
  }
}

#' @export
plot.CFscore <- function(x, ...) {
  # this plotting function should ideally be more customizable,
  # i.e. show/hide legend, colors, xlim, ylim, ....

  models <- names(x$predictions)

  plot(1, type = "n",
       xlim = c(0, 1), ylim = c(0, 1),
       asp = 1,
       xlab = "Predicted", ylab = "Observed")
  title(
    main = "Counterfactual calibration plot",
    sub = paste0("Calibration plot had everyone followed treatment option ", x$treatment_of_interest),
    col.sub = "#404040",
    cex.sub = 0.8
  )
  graphics::abline(0, 1, col = "black")
  colors <- adjustcolor(
    rep(palette()[-1], length.out = length(models)),
    alpha.f = 0.8
  )

  for (i in seq_along(models)) {
    lines(
      x = x$score$calplot[["pred", models[i]]],
      y = x$score$calplot[["obs", models[i]]],
      type = "o",
      col = colors[i],
      lw = 2
    )
  }
  legend("topleft",
         legend = models,
         col    = colors,
         lty    = 1,
         lwd    = 1,
         pch    = 1,
         bty    = "n")
  if (x$bootstrap_iterations > 0) {
    for (m in models) {
      plot(1, type = "n",
           xlim = c(0, 1), ylim = c(0, 1),
           xlab = "Predicted", ylab = "Observed",
           asp = 1)
      title(
        main = paste0("Counterfactual calibration plot for ", m),
        sub = paste0("Calibration plot had everyone followed treatment option ", x$treatment_of_interest),
        col.sub = "#404040",
        cex.sub = 0.8
      )
      for (i in 1:x$bootstrap_iterations) {
        lines(
          x = x$bootstrap$raw$calplot[[m]][[i]]$pred,
          y = x$bootstrap$raw$calplot[[m]][[i]]$obs,
          type = "o",
          col = "darkgrey"
        )
      }
      graphics::abline(0, 1, col = "black")
      lines(
        x = x$score$calplot[["pred", m]],
        y = x$score$calplot[["obs", m]],
        type = "o",
        col = "blue", lw = 2,
      )
      legend("topleft",
             legend = c("bootstrap iteration", "original (CF) data"),
             col    = c("darkgrey", "blue"),
             lty    = 1,
             lwd    = c(1,2),
             pch    = c(1,1),
             bty    = "n")
    }
  }

}


assumptions <- function(x) {

  # if ("confounders" %in% names(x)) {
  #   adjustment_text <- paste0(
  #     "{", paste0(x$ipt$confounders, collapse = ", "), "}"
  #   )
  # } else {
  #   adjustment_text <- "given IP-weights"
  # }

  pp("Estimation of the performance of the prediction model in a counterfactual
     (CF) dataset where everyone's treatment ", x$treatment_column,
     " was set to ", x$treatment_of_interest, ".")

  pp("The following assumptions must be satisfied for correct inference:")

  pp(" - Conditional exchangeability requires that the inverse probability of
  treatment weights are sufficient to adjust for confounding and selection bias
  between treatment and outcome.")

  pp("- Positivity (assess $ipt$weights for outliers)")

  pp("- Consistency (including no interference)")

  pp("- Correctly specified propensity model. Estimated treatment model is ",
     print_model(x$ipt$model), ". See also $ipt$model")

  if (x$outcome_type == "survival") {
    pp("- Censoring is accounted for with ... ")
  }
}


print_model <- function(model) {
  link <- model$family$link
  lhs_var <- model$formula[[2]]

  lhs <- paste0(link, "(", lhs_var, ")")

  var_names <- names(model$coefficients)

  coef <- round(unname(model$coefficients), 2)

  rhs <- paste(coef, var_names, sep = "*", collapse = " + ")
  rhs <- gsub("*(Intercept)", "", rhs, fixed = TRUE)
  rhs <- gsub("+ -", "- ", rhs, fixed = TRUE)

  formula <- paste(lhs, rhs, sep = " = ")
  formula
}
