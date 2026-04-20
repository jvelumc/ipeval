# Performance in observed dataset

This function computes the performance of the predictions in the given
data, which may contain a mix of treated and untreated subjects. It
exists only to demonstrate the difference between 'normal' performance
and counterfactual performance. It is not user friendly and should not
be relied on. It does not support time-to-event data.

## Usage

``` r
observed_score(
  object,
  data,
  outcome,
  metrics = c("auc", "brier", "oeratio", "calplot")
)
```

## Arguments

- object:

  One of the following three options to be validated:

  - a numeric vector, corresponding to risk predictions

  - a glm model

  - a (named) list, with one or more of the previous 2 options, for
    validating and comparing multiple models at once.

- data:

  A data.frame containing the observed outcome.

- outcome:

  The outcome, to be evaluated within data. This should be the name of a
  numeric column in data.

- metrics:

  A character vector specifying which performance metrics to be
  computed. Options are c("auc", "brier", "oeratio", "calplot").

## Value

Performance metrics in the observed dataset.
