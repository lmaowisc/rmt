# Create a recurrent event object

Create a recurrent event object

## Usage

``` r
rec(id, time, status)
```

## Arguments

- id:

  A vector of id variable.

- time:

  A vector of follow-up times.

- status:

  A vector of event type, 1 = recurrent event, 2 = death, and 0 =
  censoring;

## Value

An object of class `rec` used as an argument for
[`rmtfit`](https://lmaowisc.github.io/rmt/reference/rmtfit.md).
