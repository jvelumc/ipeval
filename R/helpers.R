is.formula <- function(x) {
  inherits(x, "formula")
}

extract_lhs <- function(data, formula) {
  mf <- model.frame(formula, data)
  unname(model.response(mf))
}

rhs_is_one <- function(formula) {
  identical(formula[[3]], 1)
}


check_missing <- function(arg) {
  is_missing <- eval(call("missing", deparse(substitute(arg))),
                     envir = parent.frame())

  if (is_missing) {
    stop("Argument ", as.name(substitute(arg)), " is missing.", call. = FALSE)
  }
}


check_missing_xor <- function(arg1, arg2) {
  is_missing1 <- eval(call("missing", deparse(substitute(arg1))),
                      envir = parent.frame())
  is_missing2 <- eval(call("missing", deparse(substitute(arg2))),
                      envir = parent.frame())

  if (!xor(is_missing1, is_missing2)) {
    stop(
      "One of arguments ", as.name(substitute(arg1)),
      " and ", as.name(substitute(arg2)),
      " must be given",
      call. = FALSE
    )
  }
}

check_input <- function(arg, class, mandatory = TRUE) {
  # check if arg is of type class, and
  # if mandatory = TRUE check if its not missing
}



simulate_time_to_event <- function(n, constant_baseline_haz, LP) {
  u <- runif(n)
  -log(u) / (constant_baseline_haz * exp(LP))
}
