# Estimate restricted mean times in favor of treatment

Estimate and make inference on the overall and component-wise restricted
mean times in favor of treatment.

## Usage

``` r
rmtfit(...)

# Default S3 method
rmtfit(id, time, status, trt, type = "multistate", ...)

# S3 method for class 'formula'
rmtfit(formula, data, ...)
```

## Arguments

- ...:

  Further arguments.

- id:

  A vector of id variable.

- time:

  A vector of follow-up times.

- status:

  For `type="multistate"`, k = entering into state \\k\\ (\\K+1\\
  represents death) and 0 = censoring; For `type="recurrent"`, 1 =
  recurrent event, 2 = death, and 0 = censoring;

- trt:

  A vector of binary variable for treatment group.

- type:

  `"multistate"` = multistate data; `"recurrent"` = recurrent event
  data.

- formula:

  A formula object. For multistate data, use `ms(id,time,status)~trt`;
  for recurrent event data, use `rec(id,time,status)~trt`.

- data:

  A data frame, which contains the variables names in the formula.

## Value

An object of class `rmtfit`. See
[`rmtfit.object`](https://lmaowisc.github.io/rmt/reference/rmtfit.object.md)
for details.

## Methods (by class)

- `rmtfit(default)`: Default

- `rmtfit(formula)`: Formula

## See also

[`rmtfit.object`](https://lmaowisc.github.io/rmt/reference/rmtfit.object.md),
[`summary.rmtfit`](https://lmaowisc.github.io/rmt/reference/summary.rmtfit.md),
[`plot.rmtfit`](https://lmaowisc.github.io/rmt/reference/plot.rmtfit.md),
[`bouquet`](https://lmaowisc.github.io/rmt/reference/bouquet.md).

## Examples

``` r
#######################
# Multistate outcome  #
#######################
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
# print the event numbers by group
obj
#> Call:
#> rmtfit.formula(formula = ms(id, time, status) ~ rx, data = colon_lev)
#> 
#>           N State 1 Death Med follow-up time
#> Control 315     177   168           5.081451
#> Lev+5FU 304     119   123           5.749487
# summarize the inference results for tau=7.5 years
summary(obj,tau=7.5)
#> Call:
#> rmtfit.formula(formula = ms(id, time, status) ~ rx, data = colon_lev)
#> 
#> Restricted mean winning time by tau = 7.5:
#>           State 1 Survival  Overall
#> Control 0.2681406 1.127625 1.395766
#> Lev+5FU 0.6140686 1.749020 2.363088
#> 
#> Restricted mean time in favor of group "Lev+5FU" by time tau = 7.5:
#>          Estimate  Std.Err Z value  Pr(>|z|)    
#> State 1  0.345928 0.072333  4.7825 1.732e-06 ***
#> Survival 0.621394 0.214220  2.9007 0.0037230 ** 
#> Overall  0.967322 0.253330  3.8184 0.0001343 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

############################
# Recurrent event outcome  #
############################
# load the HF-ACTION trial data
library(rmt)
head(hfaction)
#>        patid       time status trt_ab age60
#> 1 HFACT00001 0.60506502      1      0     1
#> 2 HFACT00001 1.04859685      0      0     1
#> 3 HFACT00002 0.06297057      1      0     1
#> 4 HFACT00002 0.35865845      1      0     1
#> 5 HFACT00002 0.39698836      1      0     1
#> 6 HFACT00002 3.83299110      0      0     1
# fit the data
obj=rmtfit(rec(patid,time,status)~trt_ab,data=hfaction)
# print the event numbers by group
obj
#> Call:
#> rmtfit.formula(formula = rec(patid, time, status) ~ trt_ab, data = hfaction)
#> 
#>     N Event 1 Event 2 Event 3 Event 4 Event 5 Event 6 Event 7 Event 8 Event 9
#> 0 221     170     117      86      56      33      23      15      13      13
#> 1 205     145      89      55      43      32      21      15      11       7
#>   Event 10 Event 11 Event 12 Event 13 Event 14 Event 15 Event 16 Event 17
#> 0       11        7        6        6        5        3        2        2
#> 1        5        4        3        2        2        2        2        2
#>   Event 18 Event 19 Event 20 Event 21 Event 22 Event 23 Event 24 Event 25
#> 0        2        1        0        0        0        0        0        0
#> 1        2        2        1        1        1        1        1        1
#>   Event 26 Death Med follow-up time
#> 0        0    57           2.390144
#> 1        1    36           2.302533
# summarize the inference results for tau=3.5 years
summary(obj,tau=3.5,Kmax=4) # aggregating results for recurrent-event
#> Call:
#> rmtfit.formula(formula = rec(patid, time, status) ~ trt_ab, data = hfaction)
#> 
#> Restricted mean winning time by tau = 3.5:
#>     Event 1   Event 2    Event 3    Event 4    Event 5    Event 6    Event 7
#> 0 0.2459671 0.1797023 0.07391981 0.05700905 0.06761022 0.03241526 0.02891853
#> 1 0.2608496 0.2245359 0.18901342 0.07904635 0.05017314 0.04026294 0.01440990
#>      Event 8    Event 9    Event 10   Event 11    Event 12     Event 13
#> 0 0.02460219 0.01354735 0.007854129 0.00103309 0.008211028 0.0003654393
#> 1 0.01168169 0.01243527 0.009208109 0.00450852 0.001553607 0.0025077249
#>       Event 14     Event 15     Event 16     Event 17     Event 18     Event 19
#> 0 0.0008724907 0.0002284718 0.0003137589 0.0001010352 0.0001884826 4.098626e-04
#> 1 0.0077205378 0.0020117607 0.0007229156 0.0003285298 0.0003447794 6.936066e-05
#>       Event 20    Event 21     Event 22     Event 23    Event 24     Event 25
#> 0 0.0001954576 0.000753361 0.0002794786 0.0005947739 0.001216739 0.0003538805
#> 1 0.0000000000 0.000000000 0.0000000000 0.0000000000 0.000000000 0.0000000000
#>     Event 26  Survival  Overall
#> 0 0.00417291 0.2960596 1.046896
#> 1 0.00000000 0.4943879 1.405772
#> 
#> Restricted mean time in favor of group "1" by time tau = 3.5:
#>           Estimate   Std.Err Z value Pr(>|z|)   
#> Event 1   0.014882  0.047748  0.3117 0.755277   
#> Event 2   0.044834  0.045834  0.9782 0.327987   
#> Event 3   0.115094  0.036098  3.1883 0.001431 **
#> Event 4+ -0.014262  0.049464 -0.2883 0.773098   
#> Survival  0.198328  0.093375  2.1240 0.033670 * 
#> Overall   0.358876  0.154388  2.3245 0.020098 * 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# frequency >=4.

```
