# A dataset from a landmark colon cancer trial

A landmark colon cancer trial on the efficacy of levamisole and
fluorouracil was reported by Moertel et al. (1990). The trial recruited
929 patients with stage C disease and randomly assigned them to
levamisole treatment alone, levamisole combined with fluorouracil, and
the control. The dataset here is restricted to the comparison between
the combined treatment and control groups, consisting of 304 and 314
patients, respectively.

## Usage

``` r
colon_lev
```

## Format

A data frame with 915 rows and 6 variables:

- id:

  Unique patient ID.

- time:

  Event time (years).

- status:

  Event type; 1 = cancer relapse, 2 = death.

- rx:

  "Lev+5FU" = combined treatment, "Control" = control.

- age:

  Patient age (years) at randomization.

- sex:

  0 = female, 1 = male.

## References

MOERTEL, C. G., FLEMING, T. R., MACDONALD, J. S., HALLER, D. G., LAURIE,
J. A., GOODMAN, P. J., UNGERLEIDER, J. S., EMERSON, W. A., TORMEY, D.
C., GLICK, J. H. et al. (1990). Levamisole and fluorouracil for adjuvant
therapy of resected colon carcinoma. New Engl. J. Med. 322, 352–358.
