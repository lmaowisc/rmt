# Estimated restricted mean times in favor of treatment

This class of objects is returned by the `rmtfit` class of functions.
Objects of this class have methods for the functions `print`, `summary`,
`plot`, and `bouquet`.

## Value

- t:

  A vector of follow-up times \\\tau\\.

- mu:

  A matrix with \\K+2\\ rows; The \\k\\th row \\(k=1,\ldots, K)\\ is
  \\\mu_k(\tau)\\, the restricted mean time in favor of treatment on the
  \\k\\th state (or recurrent event); The \\(K+1)\\th row is the net
  restricted mean survival time; The last row is the overal effect
  \\\mu(\tau)\\.

- var:

  A matrix with \\K+2\\ rows containing the variance estimates for `mu`.

- mu10, mu01:

  Matrices with \\K+1\\ rows; The \\k\\th row \\(k=1,\ldots, K)\\ is the
  restricted mean win (`mu10`) and loss (`mu01`) times by the treatment
  on the \\k\\th state (or recurrent event); The \\(K+1)\\th row is that
  on the survival time.

- ...:

## See also

[`rmtfit`](https://lmaowisc.github.io/rmt/reference/rmtfit.md),
[`summary.rmtfit`](https://lmaowisc.github.io/rmt/reference/summary.rmtfit.md),
[`plot.rmtfit`](https://lmaowisc.github.io/rmt/reference/plot.rmtfit.md),
[`bouquet`](https://lmaowisc.github.io/rmt/reference/bouquet.md).
