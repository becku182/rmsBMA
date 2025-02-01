
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rmsBMA

<!-- badges: start -->
<!-- badges: end -->

Bayesian model averaging (BMA) in the circumsence of high number of
regressors and low number of observations. It also provides standard BMA
tools. The package allows estimation of BMA statistics, performing
Extreme Bound Analysis, Bayesian model selection, jointness analysis,
and provides the user with graphical functions for model space and
coefficients.

## Installation

You can install the development version of rmsBMA from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("becku182/rmsBMA")
```

## Example 1

Example of finding robust determinants of the dependent variable in the
context of linear regression using Bayesian model averaging with reduced
model space. The example uses data on the determinants of international
trade for 11 Eurozone countries.

``` r
# loading the package
library(rmsBMA)

# inspecting the trade data set

round(head(Trade_data)[,1:10],2)

  LNTRADE LNDGEO LNRGDPPROD B L RGDPpcDIFF  GOV HUMAN   FDI  KSI
   22.45   7.01      25.16 0 1       0.09 0.03  0.23 13.74 0.29
   20.97   7.47      24.51 0 0       0.06 0.04  0.09  1.86 0.30
   23.04   7.12      26.92 0 0       0.12 0.05  0.14  1.53 0.31
   25.33   6.52      27.27 1 1       0.04 0.01  0.41  1.93 0.30
   20.62   7.64      24.55 0 0       0.21 0.03  0.40 11.81 0.50
   23.78   7.03      26.86 1 0       0.12 0.01  0.02  1.57 0.19

# First we build Post object using Super_Posterior function. We consider all the 
# models that have 6 or less regressors.

Post<-Super_Posterior(Trade_data,M=6)

# Now we can see the results obtained using uniform model prior.

round(Post[[1]],3)

              PIP     PM   PSD  PM/PSD  P(+)
 Constant   1.000  8.649 0.907   9.531 1.000
 LNDGEO     1.000 -1.222 0.089 -13.775 1.000
 LNRGDPPROD 1.000  0.880 0.029  30.183 1.000
 B          0.181  0.048 0.096   0.505 0.164
 L          0.114  0.023 0.063   0.370 0.095
 RGDPpcDIFF 0.087 -0.009 0.086  -0.102 0.965
 GOV        0.084 -0.159 0.536  -0.297 0.977
 HUMAN      0.165  0.090 0.203   0.444 0.148
 FDI        0.258  0.007 0.012   0.582 0.233
 KSI        0.364  0.876 1.115   0.786 0.344
 BCIDIFF    0.074  0.000 0.001  -0.116 0.967
 CPW        0.076  0.000 0.000   0.162 0.046
 INFVAR     0.371 -0.178 0.219  -0.816 0.981
 ARABLE     0.181  0.002 0.005   0.500 0.165
 LAND       0.234  0.000 0.000  -0.625 0.986
 ARABLEpw   0.624 -0.008 0.005  -1.437 0.990
 LANDpc     0.182  0.000 0.000   0.480 0.166
 EPCpc      0.082  0.000 0.000  -0.036 0.961

# We can compare them with the results obtained using binomial-beta model prior.

round(Post[[2]],3)

              PIP     PM   PSD  PM/PSD  P(+)
 Constant   1.000  8.463 0.837  10.113 1.000
 LNDGEO     1.000 -1.187 0.062 -19.230 1.000
 LNRGDPPROD 1.000  0.885 0.025  34.762 1.000
 B          0.080  0.020 0.066   0.302 0.072
 L          0.057  0.012 0.046   0.250 0.047
 RGDPpcDIFF 0.047 -0.009 0.063  -0.137 0.985
 GOV        0.044 -0.084 0.397  -0.211 0.988
 HUMAN      0.091  0.049 0.157   0.312 0.082
 FDI        0.089  0.002 0.008   0.254 0.072
 KSI        0.113  0.222 0.721   0.308 0.100
 BCIDIFF    0.039  0.000 0.001  -0.139 0.985
 CPW        0.038  0.000 0.000   0.107 0.022
 INFVAR     0.114 -0.046 0.142  -0.321 0.989
 ARABLE     0.101  0.001 0.004   0.347 0.092
 LAND       0.132  0.000 0.000  -0.417 0.993
 ARABLEpw   0.273 -0.003 0.005  -0.596 0.994
 LANDpc     0.102  0.000 0.000   0.337 0.093
 EPCpc      0.040  0.000 0.000  -0.055 0.982

# We might compare those results with the one obtained using Extreme Bounds Analysis.

round(Post[[2]],3)

            Lower bound Minimum   Mean Maximum Upper bound  %(+)
 Constant        -8.131  -8.131 19.445  38.905      39.566 0.880
 LNDGEO          -2.630  -2.053 -1.201   0.007       0.842 0.001
 LNRGDPPROD       0.473   0.635  0.898   1.179       1.383 1.000
 B               -0.277   0.153  1.508   2.795       4.039 1.000
 L               -3.248  -1.835  0.495   2.105       3.282 0.722
 RGDPpcDIFF      -7.607  -6.081 -2.834   0.702       2.398 0.023
 GOV            -52.379 -30.599 -3.481  11.954      31.367 0.401
 HUMAN           -5.914  -3.732 -0.146   3.741       5.876 0.440
 FDI             -0.182  -0.146 -0.035   0.128       0.188 0.296
 KSI            -17.121 -11.193 -4.077   3.542       8.761 0.135
 BCIDIFF         -0.094  -0.061 -0.003   0.087       0.147 0.442
 CPW              0.000   0.000  0.000   0.000       0.000 0.097
 INFVAR          -2.463  -1.372 -0.229   1.323       2.481 0.260
 ARABLE          -0.132  -0.070  0.000   0.072       0.148 0.531
 LAND             0.000   0.000  0.000   0.000       0.000 0.730
 ARABLEpw        -0.074  -0.050 -0.016   0.030       0.067 0.216
 LANDpc           0.000   0.000  0.000   0.000       0.000 0.104
 EPCpc            0.000   0.000  0.000   0.000       0.000 0.042

# Using Post objest we can calculate jointness measures for the regressors.

Jointness(Post)[1:5,1:5]

            LNDGEO LNRGDPPROD      B      L RGDPpcDIFF
 LNDGEO      1.000      1.000 -0.639 -0.772     -0.826
 LNRGDPPROD  1.000      1.000 -0.639 -0.772     -0.826
 B          -0.839     -0.839  1.000  0.416      0.489
 L          -0.885     -0.885  0.727  1.000      0.612
 RGDPpcDIFF -0.906     -0.906  0.753  0.796      1.000

# Next we can examine which models are the best according to PMP calculated 
# using unform model prior.

Best<-BestModels(Post,best=4)

# We can look at which variables are present in the best models 
# (PMPs are in the bottom row).

Best[[1]]

#> |           | 'No. 1' | 'No. 2' | 'No. 3' | 'No. 4' |
#> |:----------|:-------:|:-------:|:-------:|:-------:|
#> |LNDGEO     |  1.000  |  1.000  |  1.000  |  1.00   |
#> |LNRGDPPROD |  1.000  |  1.000  |  1.000  |  1.00   |
#> |B          |  0.000  |  0.000  |  0.000  |  1.00   |
#> |L          |  0.000  |  0.000  |  0.000  |  0.00   |
#> |RGDPpcDIFF |  0.000  |  0.000  |  0.000  |  0.00   |
#> |GOV        |  0.000  |  0.000  |  0.000  |  0.00   |
#> |HUMAN      |  0.000  |  0.000  |  0.000  |  0.00   |
#> |FDI        |  1.000  |  0.000  |  0.000  |  0.00   |
#> |KSI        |  1.000  |  1.000  |  1.000  |  1.00   |
#> |BCIDIFF    |  0.000  |  0.000  |  0.000  |  0.00   |
#> |CPW        |  0.000  |  0.000  |  0.000  |  0.00   |
#> |INFVAR     |  1.000  |  1.000  |  1.000  |  1.00   |
#> |ARABLE     |  0.000  |  0.000  |  0.000  |  0.00   |
#> |LAND       |  0.000  |  0.000  |  1.000  |  0.00   |
#> |ARABLEpw   |  1.000  |  1.000  |  1.000  |  1.00   |
#> |LANDpc     |  0.000  |  0.000  |  0.000  |  0.00   |
#> |EPCpc      |  0.000  |  0.000  |  0.000  |  0.00   |
#> |PMP        |  0.059  |  0.049  |  0.038  |  0.03   |

# We can also see at the best models themselves.

Best[[2]]


#> |   Model    |     'No. 1'     |     'No. 2'     |    'No. 3'     |     'No. 4'     |
#> |:----------:|:---------------:|:---------------:|:--------------:|:---------------:|
#> |  Constant  | 7.88 (1.68)***  |  9.56 (1.5)***  | 8.88 (1.52)*** | 9.18 (1.49)***  |
#> |   LNDGEO   | -1.26 (0.11)*** | -1.35 (0.1)***  | -1.36 (0.1)*** | -1.24 (0.12)*** |
#> | LNRGDPPROD |  0.9 (0.05)***  | 0.86 (0.05)***  | 0.89 (0.05)*** | 0.84 (0.05)***  |
#> |     B      |                 |                 |                |   0.29 (0.18)   |
#> |     L      |                 |                 |                |                 |
#> | RGDPpcDIFF |                 |                 |                |                 |
#> |    GOV     |                 |                 |                |                 |
#> |   HUMAN    |                 |                 |                |                 |
#> |    FDI     |  0.03 (0.02)*   |                 |                |                 |
#> |    KSI     |   2.69 (1)***   | 3.03 (1.01)***  |  3.25 (1)***   |   3.2 (1)***    |
#> |  BCIDIFF   |                 |                 |                |                 |
#> |    CPW     |                 |                 |                |                 |
#> |   INFVAR   | -0.67 (0.21)*** | -0.57 (0.21)*** | -0.6 (0.2)***  | -0.56 (0.2)***  |
#> |   ARABLE   |                 |                 |                |                 |
#> |    LAND    |                 |                 |     0 (0)*     |                 |
#> |  ARABLEpw  | -0.02 (0.01)*** |  -0.02 (0)***   |  -0.01 (0)***  |  -0.02 (0)***   |
#> |   LANDpc   |                 |                 |                |                 |
#> |   EPCpc    |                 |                 |                |                 |
#> |    R^2     |      0.96       |      0.96       |      0.96      |      0.96       |

# rmsBMA package allows drawing graphs of prior and posterior model probabilities.

# Firsly, it can be done over the models sizes.
```

<img src="man/figures/modelSizes.png" width="100%" />

``` r
# Secondly, it can be done for the best models.
```

<img src="man/figures/modelPMPs.png" width="100%" />

``` r
# We can also use rmsBMA package to depict distribiution of coeffciients over 
# the entire models space (on the condition of the inclusion 
# of the variable in a model).

# First, we can use the histogram of the coeffcients on the distance between two countries.

Coef<-coefHist(Post)
Coef$LNDGEO
```

<img src="man/figures/LnGdist.png" width="100%" />

``` r

# Secondly, we can use the kernel density of the coeffcients on 
# the GDP product of two countries.

Coef2<-coefHist(Post,kernel=1)

Coef2$LNRGDPPRO
```

<img src="man/figures/LnGDPprod.png" width="100%" />

## Example 2

One of the issues we can encounter is associated with marginal
likelihood being too close to zero in the case where we use dataset with
high number of observations. We consider how to use rmsBMA package in
this example. Again, we are finding robust determinants of the dependent
variable in the context of linear regression using Bayesian model
averaging with reduced model space. However, this time example uses data
on the determinants of international trade for 26 European Union
countries. In this context the data set consist of 325 pairs of trading
countries.

``` r
# We can try to use Super_Posterior function.

Post<-Super_Posterior(Trade_data_325)

# However we will encouter following error:

Error in Super_Posterior(Trade_data_325) : 
There is a problem with marginal likelihoods - use Super_Posterior2

#In the case when values of marginal likelihood are zero we apply Super_Posterior2 
#(You can of course start with this function).

# We haven't specified M, and in this canse the highest considered model size is equal 
# to K - total number of regressors.
# Warning: Here with 131072 models this might take a while.

Post<-Super_Posterior2(Trade_data_325)

# Again we can have a look at posterior statistics for uniform model prior.

round(Post[[1]],3)

              PIP     PM   PSD  PM/PSD  P(+)
 Constant   1.000  8.104 0.337  24.025 1.000
 B          0.847  0.361 0.067   5.417 0.846
 LNDGEO     1.000 -1.174 0.047 -24.817 1.000
 L          0.180  0.060 0.121   0.494 0.170
 LNRGDPPROD 1.000  0.889 0.009  98.771 1.000
 RGDPpcDIFF 0.072 -0.006 0.045  -0.131 0.974
 GOV        0.988 -4.042 0.560  -7.215 1.000
 HUMAN      0.056 -0.004 0.017  -0.227 0.979
 CPW        0.157  0.000 0.000   0.453 0.146
 INFVAR     0.941 -0.038 0.005  -6.943 1.000
 ARABLE     0.229  0.001 0.002   0.610 0.221
 ARABLEpw   0.754 -0.007 0.002  -3.072 0.997
 LAND       0.061  0.000 0.000   0.119 0.036
 LANDpc     0.753  0.000 0.000   3.217 0.750
 EPCpc      0.294  0.000 0.000   0.706 0.283
 FDI        0.868  0.025 0.007   3.629 0.867
 KSI        0.423 -0.432 0.408  -1.060 0.994
 BCIDIFF    0.063  0.000 0.001  -0.123 0.976

# We might compare those results with the one obtained using Extreme Bounds Analysis.

round(Post[[3]],3)

            Lower bound Minimum   Mean Maximum Upper bound  %(+)
 Constant        -3.250  -3.250 15.295  33.176      33.176 0.871
 B               -0.070   0.235  1.220   2.728       3.450 1.000
 LNDGEO          -1.882  -1.568 -1.185  -0.745      -0.355 0.000
 L               -1.002   0.066  1.042   3.285       4.388 1.000
 LNRGDPPROD       0.795   0.847  0.897   0.995       1.065 1.000
 RGDPpcDIFF      -4.644  -3.616 -0.709   0.588       1.501 0.118
 GOV            -26.664 -20.093 -6.677   0.223       3.981 0.000
 HUMAN           -2.380  -1.137 -0.209   0.888       2.011 0.181
 CPW              0.000   0.000  0.000   0.000       0.000 0.710
 INFVAR          -0.226  -0.158 -0.045   0.033       0.106 0.036
 ARABLE          -0.051  -0.030 -0.003   0.015       0.025 0.506
 ARABLEpw        -0.093  -0.076 -0.028   0.014       0.022 0.112
 LAND             0.000   0.000  0.000   0.000       0.000 0.605
 LANDpc           0.000   0.000  0.000   0.000       0.000 0.319
 EPCpc            0.000   0.000  0.000   0.000       0.000 0.820
 FDI             -0.160  -0.120  0.023   0.135       0.181 0.762
 KSI            -12.025 -10.027 -3.469   0.359       2.640 0.001
 BCIDIFF         -0.062  -0.037  0.013   0.061       0.086 0.695

# In the case of M=K, we have really informative graphs depicting prior and 
# posterior model probabilities over the model sizes.

modelSizes(Post)
```

<img src="man/figures/ModelSizes325.png" width="100%" />

``` r
# In this case it is also very instructive to have a look a prior and posterior 
# probabilities of the best 200 models.

modelPMP(Post,Top=200)
```

<img src="man/figures/ModelPMP325.png" width="100%" />

``` r
# Lets now have a look a the histogram of the coefficients of the product of GDP of 
# the trading countries (LNRGDPPROD).

Coef<-coefHist(Post)
Coef$LNRGDPPROD
```

<img src="man/figures/LnPROD_325.png" width="100%" />

``` r
# Now we use the kernel density of the coeffcients on the distance 
# between two countries (LNDGEO).

Coef2<-coefHist(Post,kernel=1)
Coef2$LNDGEO
```

<img src="man/figures/LnGDIST_325.png" width="100%" />

``` r
# Finally, lets have a look at the 4 best models according to PMPs based 
# on uniform model prior.

Best<-BestModels(Post)
Best[[2]]

#> |   Model    |     'No. 1'     |     'No. 2'     |     'No. 3'     |     'No. 4'     |
#> |:----------:|:---------------:|:---------------:|:---------------:|:---------------:|
#> |  Constant  |  8.1 (0.68)***  | 8.12 (0.67)***  | 8.19 (0.67)***  | 7.82 (0.69)***  |
#> |     B      | 0.44 (0.14)***  | 0.41 (0.14)***  | 0.41 (0.14)***  | 0.46 (0.14)***  |
#> |   LNDGEO   | -1.19 (0.06)*** | -1.15 (0.07)*** | -1.1 (0.07)***  | -1.18 (0.06)*** |
#> |     L      |                 |                 |                 |                 |
#> | LNRGDPPROD | 0.89 (0.02)***  | 0.89 (0.02)***  | 0.88 (0.02)***  | 0.89 (0.02)***  |
#> | RGDPpcDIFF |                 |                 |                 |                 |
#> |    GOV     | -4.33 (0.86)*** | -3.59 (0.92)*** | -3.47 (0.93)*** | -4.35 (0.86)*** |
#> |   HUMAN    |                 |                 |                 |                 |
#> |    CPW     |                 |                 |                 |                 |
#> |   INFVAR   | -0.04 (0.01)*** | -0.03 (0.01)*** | -0.04 (0.01)*** | -0.04 (0.01)*** |
#> |   ARABLE   |                 |                 |                 |    0.01 (0)*    |
#> |  ARABLEpw  |  -0.01 (0)***   |  -0.01 (0)***   |  -0.01 (0)***   |  -0.01 (0)***   |
#> |    LAND    |                 |                 |                 |                 |
#> |   LANDpc   |    0 (0)***     |    0 (0)***     |                 |    0 (0)***     |
#> |   EPCpc    |                 |                 |    0 (0)***     |                 |
#> |    FDI     | 0.03 (0.01)***  | 0.03 (0.01)***  | 0.03 (0.01)***  | 0.03 (0.01)***  |
#> |    KSI     |                 | -0.87 (0.41)**  | -1.12 (0.42)*** |                 |
#> |  BCIDIFF   |                 |                 |                 |                 |
#> |    R^2     |      0.93       |      0.93       |      0.93       |      0.93       |
```
