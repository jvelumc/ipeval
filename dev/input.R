# library(rlang)
# library(lobstr)
library(survival)

df <- data.frame(X = runif(5), Y = rbinom(5, 1, 0.5))
model <- glm(Y ~ X, family = "binomial", data = df)

Y <- 1

myfunc <- function(data, object, outcome) {

  object <- make_list_if_not_list(object)


  object_name <- substitute(object)

  var <- eval(substitute(outcome), data)
  list(
    "names" = object_name,
    "outcome" = var,
    "object" = object
  )

}

Surv(df$X, df$Y)

myfunc(df, object = c(1,2,3), Q)


View(expr(f(x, "y", 1)))

enexpr()

ast(list(a, b, c))

local()


ast(1+2+3)
