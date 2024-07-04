# PMP.jl
*PMP.jl* is a package to ease the estimation of the probable maximum precipitation (PMP) through various methods.


[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build status](https://github.com/JuliaExtremes/PMP.jl/workflows/CI/badge.svg)](https://github.com/JuliaExtremes/PMP.jl/actions)
[![codecov](https://codecov.io/gh/JuliaExtremes/PMP.jl/branch/master/graph/badge.svg?token=d8ecbbb8-ea4f-42cc-8317-39e8ceb648fb)](https://codecov.io/gh/JuliaExtremes/PMP.jl)
[![documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaextremes.github.io/PMP.jl/stable/)
[![documentation latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliaextremes.github.io/PMP.jl/dev/)

The [World Meteorological Organization (WMO) guide](https://library.wmo.int/index.php?lvl=notice_display&id=1302#.ZLlVeezMKeA) (2009) define the PMP as "the greatest depth of precipitation for a given duration meteorologically possible for a design watershed or a given storm area at a particular location at a particular time of year, with no allowance made for long-term climatic trends".

## Documentation
### Tutorial
The documentation includes a tutorial on how to compute the PMP using different approaches:

- Moisture maximization (MoistureMaximization.md),
- Hershfield empirical method (OtherMethods.md),
- Univariate GEV method (OtherMethods.md),
- Beta method (BetaMethod.md).

The WMO (2009) describes in detail steps of the first two methods. The two other ones are statistical methods developped to alleviate flaws raised by the [scientific community](https://nap.nationalacademies.org/catalog/27460/modernizing-probable-maximum-precipitation-estimation).

### Notebooks
The documentation also includes julia notebooks that tests the beta method, which is a novel approach developped here. The *01-4params.ipynb* and *02-3params.ipynb* files test the 4 and 3 parameters Pearson type I models respectively, with standard parametrization and frequentist methods (method of moments, MLE, mixed method). The *03-3params-bayes* folder gathers samples and files used to test two bayesian methods (Gibbs sampler, NUTS), and *04-reparam.ipynb* tests the 3 parameters reparametrized Pearson type I model. Finally, *05-rain-data.ipynb* and *06-usual-PMP-methods.ipynb* compute PMP estimations on rain datasets observed at two stations (Pierre-Elliott-Trudeau airport of Montréal and St-Hubert) with the Beta method for the first file and moisture maximization and Hershfield empirical method for the second.
