# PMP estimation with other methods

This example shows how to estimate 72h PMP by univariate GEV and Hershfield methods with *PMP.jl* and the recommended methodology by the [World Meteorological Organization](https://library.wmo.int/index.php?lvl=notice_display&id=1302#.ZLlVeezMKeA) (WMO). The example uses observed data from the Montréal-Trudeau International Airport station. To avoid solid precipitation, we only consider data from May to October and from 1953 to 2012. 


## Load required Julia packages

Before executing this tutorial, make sure to have installed the following packages:

- *DataFrames.jl* (for using the DataFrame type)
- *PMP.jl*

and import them using the following command:
 ```@repl OtherMethods
using DataFrames, PMP
```


## Load required dataset

Loading the observed daily precipitations (in mm) :
```@example OtherMethods
# Load the data
rain = PMP.dataset("rain")
println("") # hide
```

---
## Univariate GEV method

The univariate GEV method is not a WMO-recommended procedure. It estimates the return value of the rain for a given return period. In this example, we choose a 100-year return period.

```@example OtherMethods
pmp_gev = PMP_GEV(rain.Rain, rain.Date, 100, 24, 72)
println("") # hide
```

We could estimate the 24h PMP with the following :
```@example OtherMethods
pmp_gev = PMP_GEV(rain.Rain, rain.Date, 100)
println("") # hide
```

---
## Hershfield method

The Hershfield method is an empirical approach using the mean and standard deviation of annual maximum precipitations and an abstracted statistic `K`. The fourth chapter of the WMO manual covers this estimation method in more detail. This present package does not take into consideration the size of the dataset. `PMP_Hershfield` function can take a chosen `K` in argument or calculate one with the dataset :

```@example OtherMethods
pmp_hershfield_k15 = PMP_Hershfield(rain.Rain, rain.Date, 15, 24, 72)
pmp_hershfield = PMP_Hershfield(rain.Rain, rain.Date, 24, 72)
println("") # hide
```

We could estimate the 24h PMP with the following :
```@example OtherMethods
pmp_hershfield_k15 = PMP_Hershfield(rain.Rain, rain.Date, 15)
pmp_hershfield = PMP_Hershfield(rain.Rain, rain.Date)
println("") # hide
```