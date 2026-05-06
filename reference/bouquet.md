# Bouquet plot

Construct the bouquet plot based on the estimated stage-wise restricted
mean win/loss times.

## Usage

``` r
bouquet(
  x,
  Kmax = NULL,
  xlim = NULL,
  ylim = NULL,
  xlab = "Restricted mean win/loss time",
  ylab = "Follow-up time",
  group.label = TRUE,
  cex.group = 1,
  ...
)
```

## Arguments

- x:

  An object returned by
  [`rmtfit`](https://lmaowisc.github.io/rmt/reference/rmtfit.md).

- Kmax:

  A positive integer; If specified, the stage-wise estimates over
  `Kmax`\\,\ldots,K\\ will be aggregated.

- xlim:

  The x limits of the plot.

- ylim:

  The y limits of the plot.

- xlab:

  A label for the x axis, defaults to a description of x.

- ylab:

  A label for the y axis, defaults to a description of x.

- group.label:

  If `TRUE`, group labels will appear on the two sides of the plot.

- cex.group:

  Font size of the group labels if `group.label=TRUE`.

- ...:

  Other arguments that can be passed to the underlying `plot` method.

## Value

No return value, called for side effects.

## See also

[`rmtfit`](https://lmaowisc.github.io/rmt/reference/rmtfit.md),
[`summary.rmtfit`](https://lmaowisc.github.io/rmt/reference/summary.rmtfit.md),
[`plot.rmtfit`](https://lmaowisc.github.io/rmt/reference/plot.rmtfit.md).

## Examples

``` r
# load the colon cancer trial data
library(rmt)
head(colon_lev)
#>   id      time status      rx sex age
#> 1  1 2.6502396      1 Lev+5FU   1  43
#> 2  1 4.1642710      2 Lev+5FU   1  43
#> 3  2 8.4517454      0 Lev+5FU   1  63
#> 4  3 1.4839151      1 Control   0  71
#> 5  3 2.6365503      2 Control   0  71
#> 6  4 0.6707734      1 Lev+5FU   0  66
# fit the data
obj=rmtfit(ms(id,time,status)~rx,data=colon_lev)
# bouquet plot
bouquet(obj)
```
