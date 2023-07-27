
# PMP estimation by moisture maximization 

This example shows how to estimate 72h PMP by moisture maximization with *PMP.jl* and the recommended methodology by the [World Meteorological Organization](https://library.wmo.int/index.php?lvl=notice_display&id=1302#.ZLlVeezMKeA) (WMO). The example uses observed data from the Montréal-Trudeau International Airport station. To avoid solid precipitation, we only consider data from May to October and from 1953 to 2012. 


## Load required Julia packages

Before executing this tutorial, make sure to have installed the following packages:

- *DataFrames.jl* (for using the DataFrame type)
- *PMP.jl*

and import them using the following command:
 ```@repl MoistureMaximization
using DataFrames, PMP
```


## Load required datasets

Loading the observed daily precipitations (in mm) and observed hourly dew point (in °C):
```@example MoistureMaximization
# Load the data
rain = PMP.dataset("rain")
dewpoint = PMP.dataset("dewpoint")
 
println("") # hide
```

---
## Storm selection

First, we need to select storms for maximization. We want to focus on maximizing the 10% largest precipitation event of each year.

```@example MoistureMaximization
# Select "PMP magnitude" storms
p = 0.1    # 10% 
d1 = 24    # One observation per day
d2 = 72    # PMP duration

storm = storm_selection(rain.Rain, rain.Date, p, d1, d2)
println("") # hide
```

We also could have used `total_precipitation` then `storm_selection` as follow :

```@example MoistureMaximization
# Aggregate precipitations of duration d1 to precipitations of duration d2 without overlap
rain_on_72h = total_precipitation(rain.Rain, rain.Date, d1, d2)

# Select the p greatest storms of each year
storm = storm_selection(rain_on_72h.Rain, rain_on_72h.Date, p)
println("") # hide
```

---
## Precipitable water calculation

At Montréal-Trudeau International Airport station, we do not have direct access to precipitable water (PW) observations; instead, we need to convert dewpoint (in °C) to PW (in mm). The dewpoint is measured hourly, and the highest persisting dewpoint represents the data that best characterizes the humidity for each storm. According to the WMO manual, the persisting dewpoint for a specified time interval is defined as 'the value equalled or exceeded at all observations during the period.' For our purposes, we will use 12-hour time intervals, as recommended by the WMO. Furthermore, we are specifically interested in dewpoint observations taken during the PMP magnitude storms determined in the previous section. Accordingly, we need to select those relevant observation from the dewpoint dataset. 

The `PW_storm` function utilizes storm dates, the dewpoint dataset, PMP duration, dewpoint observation frequency, and the selected time interval for persisting dewpoint to return the PW of each storm of interest.

```@example MoistureMaximization
obs_freq = 1     # One observation per hour
time_int = 12    # Keep the smallest observation on each moving window of size 12h

# Calculate appropriate pw for each storm
pw_storm = PW_storm(storm.Date, dewpoint.Dew, dewpoint.DateTime, d2, obs_freq, time_int)
println("") # hide
```

If the dewpoint dataset contains observations associated only with the storms of interest, it is possible to use `get_max_persisting_dew` followed by `dewpoint_to_PW` to calculate the PW for each precipitation event. Some methods ([Beauchamp *et al.*, 2013](https://doi.org/10.1002/wrcr.20336), [Rousseau *et al.*, 2014](https://doi.org/10.1016/j.jhydrol.2014.10.053)) suggest using humidity data from before the event in addition to the observations taken during the storm. Hence, it could be appropriate to either utilize these two functions or modify `PW_storm` to take those data into consideration.

---
## Maximization ratio and PMP estimation

To maximize the PMP magnitude storms, the precipitation has to be multiplied by the maximization ratio (MR), which is defined by ``\frac{PW_{max}}{PW_{storm}}``. ``PW_{max}`` can be the maximum precipitable water recorded during the month of the event of the entire dataset or the return value for a given return period for the month of the event. We calculate both fot all the months of interest (May to October): 

```@example MoistureMaximization
# PW dataset construction
pw_dataset = DataFrame(Date = dewpoint.DateTime)
pw_dataset.PW = dewpoint_to_PW.(dewpoint.Dew)

pw_max = PW_max(pw_dataset.PW, pw_dataset.Date)                  # PW max
pw_100 = PW_return_period(pw_dataset.PW, pw_dataset.Date, 100)   # PW 100

println("") # hide
```

If we want information (date, maximization ratio, effective precipitation (EP) and maximized rain) on all the maximized storms, we can call the `storm_maximization` function : 

```@example MoistureMaximization
maximized_storm_PW_max = storm_maximization(storm.Rain, pw_storm.PW, storm.Date, pw_max.PW_max)
maximized_storm_PW_100 = storm_maximization(storm.Rain, pw_storm.PW, storm.Date, pw_100.PW_rp)
println("") # hide
```

Alternatively, we could exclusively estimate the PMP as follows :

```@example MoistureMaximization
pmp_max = PMP_mm(storm.Rain, pw_storm.PW, storm.Date, pw_max.PW_max)
pmp_100 = PMP_mm(storm.Rain, pw_storm.PW, storm.Date, pw_100.PW_rp)
println("") # hide
```


