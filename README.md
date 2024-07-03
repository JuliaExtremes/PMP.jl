# PMP.jl
Probable Maximum Precipitation (PMP) estimation


[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build status](https://github.com/JuliaExtremes/PMP.jl/workflows/CI/badge.svg)](https://github.com/JuliaExtremes/PMP.jl/actions)
[![codecov](https://codecov.io/gh/JuliaExtremes/PMP.jl/branch/master/graph/badge.svg?token=d8ecbbb8-ea4f-42cc-8317-39e8ceb648fb)](https://codecov.io/gh/JuliaExtremes/PMP.jl)
[![documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaextremes.github.io/PMP.jl/stable/)
[![documentation latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliaextremes.github.io/PMP.jl/dev/)

The [World Meteorological Organization (WMO) guide](https://library.wmo.int/index.php?lvl=notice_display&id=1302#.ZLlVeezMKeA) (2009) define the PMP as "the greatest depth of precipitation for a given duration meteorologically possible for a design watershed or a given storm area at a particular location at a particular time of year, with no allowance made for long-term climatic trends".

*PMP.jl* is a package to ease the estimation of the PMP through various methods:

- Moisture maximization,
- Hershfield empirical method,
- Univariate GEV method,
- Beta method.

The WMO (2009) describes in detail steps of the first two methods. The two other ones are statistical methods developped to alleviate flaws raised by the [scientific community](https://nap.nationalacademies.org/catalog/27460/modernizing-probable-maximum-precipitation-estimation).
