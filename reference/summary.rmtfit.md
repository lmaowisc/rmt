# Summary of the analysis results

Summarize the overall and stage-wise inferential results for the
restricted mean times in favor of treatment at a user-specified length
of follow-up.

## Usage

``` r
# S3 method for class 'rmtfit'
summary(object, tau = NULL, Kmax = NULL, ...)
```

## Arguments

- object:

  An object returned by
  [`rmtfit`](https://lmaowisc.github.io/rmt/reference/rmtfit.md).

- tau:

  A positive real number for the follow-up time; Default is the maximum
  event time in the data.

- Kmax:

  A positive integer; If specified, the stage-wise estimates over
  `Kmax`\\,\ldots,K\\ will be aggregated.

- ...:

  Additional arguments affecting the summary produced.

## Value

An object of class `summary.rmtfit` with components

- WL:

  A \\2\times(K+2)\\-dimensional matrix; Each row contains the estimates
  for the stage-wise and overall restricted mean win times for each
  group.

- tab:

  A \\(K+2)\times 4\\-dimensional matrix summarizing the inferential
  results for the stage-wise and overall restricted mean times in favor
  of treatment; Columns include `Estimate`, `Std.Err`, `Z value`, and
  `Pr(>|z|)`.

- ...:

## See also

[`rmtfit`](https://lmaowisc.github.io/rmt/reference/rmtfit.md),
[`plot.rmtfit`](https://lmaowisc.github.io/rmt/reference/plot.rmtfit.md),
[`bouquet`](https://lmaowisc.github.io/rmt/reference/bouquet.md).

## Examples

``` r
#See examples for rmtfit().
```
