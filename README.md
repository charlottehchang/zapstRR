# zapstRR: ZoologicAl Package (Zap!) for Randomized Response Technique Analysis

The package `zapstRR` provides univariate (single item) and multivariate (multiple item) estimators for randomized response technique (RRT) data. We aim to improve the functionality of RRT survey design and analysis by including the following functions:

* `RRunivariate`: analysis (prevalence and logistic regression) of single RRT items
* `RRsumscore`: prevalence and ordinal regression of the sums cores of multiple RRT items
* `RRirt`: multivariate (two or more RRT items) prevalence and regression analysis
⋅⋅+ `RRirt` can also be used to estimate evasive response bias--or the proportion of individuals that consistently refuse to answer yes honestly or when prompted by the randomizer.

The three functions above take the same input call: `(y = Response Data, x = Predictor Data, p00 = Probability of observing No conditional on a true No, p11 = Probability of observing Yes conditional on a true Yes)`. Note that if nothing is provided for `x`, then the functions will return proportions (i.e. the percentage of respondents who exhibited a specific sensitive trait, or possessed some ordinal sum of traits).

Additionally, `zapstRR` contains these functions:

* `RRbootstrap`: parameter estimation using non-parametric bootstrapping for all three models: `univariate`, `sumscore`, and `irt`.
    + `RRboostrap` permits users to estimate confidence intervals when there is insufficient information to otherwise use maximum likelihood estimation to determine the standard error of a parameter.
* `powerRRT`: power analysis of study designs
* `RRsimulate`: simulation of RRT data under the three models: `univariate`, `sumscore`, and `irt`.

zapstRR has two example datasets:

* `huntRRdata`: RRT data on barbet (**Psilopogon** spp.), bulbul (**Pycnonotus**, **Hypsipetes**, **Iole**, and **Ixos** spp.), and partridge (**Arborophila** spp.) hunting with respondent age (standardized to z-scores) as a covariate (Chang **et al.** 2017, in prep.)
* `RRdata`: illicit or sensitive behaviors among university undergraduates (e.g. bullying, drug usage) (adapted from Cobo Rodriguez **et al.** 2015; please see their package [RRTCS](https://cran.r-project.org/web/packages/RRTCS/index.html) for more information). 

## Downloading

At present, the package can be downloaded using [devtools](https://cran.r-project.org/web/packages/devtools/) using the command 
```r
library("devtools")
install_github("charlottehchang/zapstRR")
``` 