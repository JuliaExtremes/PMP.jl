
# PMP estimation by moisture maximization 

This example shows how to estimate 72h PMP by moisture maximization with *PMP.jl* and the recommended methodology by the [World Meteorological Organization](https://library.wmo.int/index.php?lvl=notice_display&id=1302#.ZLlVeezMKeA). The example uses observed data of the Montréal-Trudeau station. To avoid solid precipitation, we only consider data from May to October and from 1953 to 2012. 

## Load required Julia packages

Before executing this tutorial, make sure to have installed the following packages:

- *CSV.jl* (for loading the data)
- *DataFrames.jl* (for using the DataFrame type)
- *Distributions.jl* (for using distribution objects)
- *Dates.jl*
- *Extremes.jl*
- *PMP.jl*

and import them using the following command:
 ```@repl stationary
import Pkg; Pkg.add("CSV"), Pkg.add("Extremes")
using CSV, DataFrames, Dates, Distributions
using Extremes, PMP
```

## Storm selection

First we need to select storms to be maximized. As we have observed daily precipitations (in mm) and the desired PMP is 72h, let `d_1` = 24 and `d_2` = 72.

Loading the observed daily precipitations (in mm)
```@example stationary
# Load the data 
 
println("") # hide
```

## Precipitable water calculation

## Maximization ratio and PMP estimation