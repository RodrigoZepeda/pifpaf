
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `pifpaf`: Estimation of Potential Impact Fractions (pif) and Population Attributable Fractions (paf) <img src="man/figures/logo.png"  style="float:right; height:135px;" alt="" />

**THIS IS STILL A WORK IN PROCESS**

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/last-week/pifpaf)](https://cran.r-project.org/package=pifpaf)
[![Codecov test
coverage](https://codecov.io/gh/RodrigoZepeda/pifpaf/branch/main/graph/badge.svg)](https://app.codecov.io/gh/RodrigoZepeda/pifpaf?branch=main)
[![CRAN
status](https://www.r-pkg.org/badges/version/pifpaf)](https://CRAN.R-project.org/package=pifpaf)
<!-- badges: end -->

> The `pifpaf` package corresponds to an **update** on the previous
> [homonimous package](https://github.com/INSP-RH/pifpaf) developed at
> INSP.

# Installation

You can install the development version of `pifpaf` from
[GitHub](https://github.com/pifpaf) with:

``` r
#install.packages("remotes")
remotes::install_github("RodrigoZepeda/pifpaf", dependencies = TRUE)
```

# Usage

The main purpose of the package is to calculate **potential impact
fractions** and **population attributable fractions** for a specific
case: when the exposure information comes from cross-sectional studies
(in particular, surveys) and the relative risk information from the
literature (*e.g.* metanalysis).

In our package we assume the user has access to the data either as
**individual-level** exposure (*e.g.* the user has access to the survey)
or the user has access to **aggregated data** on the exposure (*e.g.*
the user has the population-mean exposure and its confidence interval).
The following table shows which function to use depending on the case:

| **Exposure data**                          | Population Attributable Fraction                                   | Potential Impact Fraction                                   |
|--------------------------------------------|--------------------------------------------------------------------|-------------------------------------------------------------|
| [Individual-level](#individual-level-data) | [`paf`](#population-attributable-fraction)                         | [`pif`](#potential-impact-fraction)                         |
| [Aggregated](#aggregated-data)             | [`paf_approximate`](#approximate-population-attributable-fraction) | [`pif_approximate`](#approximate-potential-impact-fraction) |

> **Note** If the user has available both **individual-level** and
> **aggregated data** they should prefer the **individual-level**
> information as it captures better the sample’s variability. The
> population attributable fraction and potential impact fractions
> estimated from **aggregated data** are only approximations and are
> biased estimates.

## Individual-level data

In the case of individual-level data the Population Attributable
Fraction (PAF) and the Potential Impact Fraction (PIF) can be estimated
directly via the `pif()` and `paf()` functions. For the purpose of this
example, we will use data from the Mexican National Health and Nutrition
Survey of 2018 (ENSANUT 2018):

``` r
data(ensanut)

head(ensanut)
#>   age    sex    weight strata hypertension delta_na_phase_1
#> 1  28 Female  32133.58    222        FALSE         5.899257
#> 2  24 Female  75955.15    222        FALSE        19.112181
#> 3  45 Female  16075.99    222        FALSE        49.546444
#> 4  39   Male  83462.51    222        FALSE         6.273072
#> 5  41   Male  84734.87    223        FALSE         0.000000
#> 6  61   Male 118740.93    223        FALSE         0.000000
#>   systolic_blood_pressure age_group
#> 1                   118.0   [25,30)
#> 2                    92.0   [20,25)
#> 3                    93.0   [45,50)
#> 4                   116.0   [35,40)
#> 5                   110.5   [40,45)
#> 6                   122.5   [60,65)
```

The data contains information on the age, sex, weight (grams),
hypertension status, systolic blood pressure and change in sodium from a
policy intervention.

### Population Attributable Fraction

### Potential Impact Fraction

### Model’s diagnostics

## Aggregated data

### Approximate Population Attributable Fraction

### Approximate Potential Impact Fraction

**WORK IN PROGRESS**
