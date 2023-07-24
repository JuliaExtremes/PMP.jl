
# PMP estimation by moisture maximization 

This example shows how to estimate 72h PMP by moisture maximization with *PMP.jl* and the recommended methodology by the [World Meteorological Organization](https://library.wmo.int/index.php?lvl=notice_display&id=1302#.ZLlVeezMKeA). The example uses observed data of the Montréal-Trudeau International Airport station. To avoid solid precipitation, we only consider data from May to October and from 1953 to 2012. 

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
using CSV, DataFrames, Dates, Distributions
using Extremes, PMP
```

## Load required datasets

Loading the observed daily precipitations (in mm) and observed hourly dew point (in °C)
```@example stationary
# Load the data
rain = PMP.dataset("rain")
dewpoint = PMP.dataset("dewpoint")
 
println("") # hide
```

## Storm selection

First we need to select storms to be maximized. We only want to maximize the 10\% biggest storms of each year.

```@example stationary
# Select "PMP magnitude" storms
p = 0.1    # 10% 
d1 = 24    # frenquency of the observations
d2 = 72    # duration of PMP

storms = PMP.storm_selection(rain.Rain, rain.Date, p, d1, d2)
println("") # hide
```

We also could have used `total_precipitation` then `storm_selection` as follow :

```@example stationary
p = 0.1    # 10% 
d1 = 24    # frenquency of the observations
d2 = 72    # duration of PMP

rain_on_72h = PMP.total_precipitation(rain.Rain, rain.Date, d1, d2)
storms = PMP.storm_selection(rain_on_72h.Rain, rain_on_72h.Date, p)
println("") # hide
```

## Precipitable water calculation

At Montréal-Trudeau International Airport station, we do not have access to precipitable water observations. We need to convert dewpoint (in °C) to precipitable water (in mm). However, as the dewpoint is measured hourly, we also need to determine which data will represent humidity for each storm. In its manual, the WMO has defined the highest persisting dewpoint for some specified time interval as "the value equalled or exceeded at all observations during the period". We will here use a 12h period as it is what recommended by the WMO and our time interval is the duration of each storm (72h). 

We first need to select dewpoint data for each storms of magnitude PMP selected in the previous section
```@example stationary
storm_dewpoint = 

println("") # hide
```

## Maximization ratio and PMP estimation