# Create a multistate event object

Create a multistate event object

## Usage

``` r
ms(id, time, status)
```

## Arguments

- id:

  A vector of id variable.

- time:

  A vector of follow-up times.

- status:

  A vector of event type, `k` if transitioning to state \\k\\, 0 if
  censored and \\K+1\\ represents death.

## Value

An object of class `ms` used as an argument for
[`rmtfit`](https://lmaowisc.github.io/rmt/reference/rmtfit.md).
