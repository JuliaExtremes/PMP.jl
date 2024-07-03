# PMP estimation by beta method

This example shows how to estimate 24h PMP by beta method with *PMP.jl*. The example uses observed daily precipitation data from the Montréal-Trudeau International Airport station from 1953 to 2012. To avoid solid precipitation, we only consider data from May to October. 

## Load required Julia packages

Before executing this tutorial, make sure to have installed the following packages:

- *DataFrames.jl* (for using the DataFrame type)
- *PMP.jl*

and import them using the following command:
 ```@repl BetaMethod
using DataFrames, PMP
```


## Load required dataset

Loading the observed daily precipitations (in mm) :
```@example BetaMethod
# Load the data
rain = PMP.dataset("rain")
rain = filter(:Rain => n -> n>0, rain); # only precipitation data
println("") # hide
```

---
## Beta method

The Beta method is not a WMO-recommended procedure. It estimates the upper limit of the beta distribution fitted to rain data. The method of moment, mixed method and Gibbs sampler are recommended over the maximum likelihood estimation. The four parameters with fixed lower limit is also prefered to the other models. Complete list of possibilities for this method can be found in the index.

Using moment method:
```@example BetaMethod
pmp_moment = fit_mme(PearsonType1, rain.Rain, 0.2).b
println("") # hide
```

Using mixed method:
```@example BetaMethod
pmp_mixed = getinitialvalues(PearsonType1, rain.Rain, 0.2)[2]
println("") # hide
```

Using Gibbs sampler (returns the Markov chain):
```@example BetaMethod
pmp_gibbs = fit_bayes_MH(PearsonType1, rain.Rain, 0.2)[1]
println("") # hide
```
