var documenterSearchIndex = {"docs":
[{"location":"MoistureMaximization/#PMP-estimation-by-moisture-maximization","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"","category":"section"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"This example shows how to estimate 72h PMP by moisture maximization with PMP.jl and the recommended methodology by the World Meteorological Organization. The example uses observed data from the Montréal-Trudeau International Airport station. To avoid solid precipitation, we only consider data from May to October and from 1953 to 2012. ","category":"page"},{"location":"MoistureMaximization/#Load-required-Julia-packages","page":"PMP estimation by moisture maximization","title":"Load required Julia packages","text":"","category":"section"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"Before executing this tutorial, make sure to have installed the following packages:","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"CSV.jl (for loading the data)\nDataFrames.jl (for using the DataFrame type)\nDistributions.jl (for using distribution objects)\nDates.jl\nExtremes.jl\nPMP.jl","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"and import them using the following command:","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"using CSV, DataFrames, Dates, Distributions\nusing Extremes, PMP","category":"page"},{"location":"MoistureMaximization/#Load-required-datasets","page":"PMP estimation by moisture maximization","title":"Load required datasets","text":"","category":"section"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"Loading the observed daily precipitations (in mm) and observed hourly dew point (in °C)","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"# Load the data\nrain = PMP.dataset(\"rain\")\ndewpoint = PMP.dataset(\"dewpoint\")\n \nprintln(\"\") # hide","category":"page"},{"location":"MoistureMaximization/#Storm-selection","page":"PMP estimation by moisture maximization","title":"Storm selection","text":"","category":"section"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"First, we need to select storms for maximization. We want to focus on maximizing the 10% largest precipitation event of each year.","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"# Select \"PMP magnitude\" storms\np = 0.1    # 10% \nd1 = 24    # frenquency of the observations\nd2 = 72    # duration of PMP\n\nstorm = PMP.storm_selection(rain.Rain, rain.Date, p, d1, d2)\nprintln(\"\") # hide","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"We also could have used total_precipitation and then storm_selection as follow :","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"rain_on_72h = PMP.total_precipitation(rain.Rain, rain.Date, d1, d2)\nstorm = PMP.storm_selection(rain_on_72h.Rain, rain_on_72h.Date, p)\nprintln(\"\") # hide","category":"page"},{"location":"MoistureMaximization/#Precipitable-water-calculation","page":"PMP estimation by moisture maximization","title":"Precipitable water calculation","text":"","category":"section"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"At Montréal-Trudeau International Airport station, we do not have direct access to precipitable water (PW) observations; instead, we need to convert dewpoint (in °C) to PW (in mm). The dewpoint is measured hourly, and the highest persisting dewpoint represents the data that best characterizes the humidity for each storm. According to the WMO manual, the persisting dewpoint for a specified time interval is defined as 'the value equalled or exceeded at all observations during the period.' For our purposes, we will use 12-hour time intervals, as recommended by the WMO. Furthermore, we are specifically interested in dewpoint observations taken during the PMP magnitude storms determined in the previous section. Accordingly, we need to select those relevant observation from the dewpoint dataset. ","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"The storm_PW function utilizes storm dates, the dewpoint dataset, PMP duration, dewpoint observation frequency, and the selected time interval for persisting dewpoint to return the PW of each storm of interest.","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"pw_storm = storm_PW(storm.Date, dewpoint.Dew, dewpoint.DateTime, 72, 1, 12)\nprintln(\"\") # hide","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"If the dewpoint dataset contains observations associated only with the storms of interest, it is possible to use get_max_persisting_dew followed by dewpoint_to_PW to calculate the PW for each precipitation event. Some methods (Beauchamp et al., 2013, Rousseau et al., 2014) suggest using humidity data from before the event in addition to the observations taken during the storm. Hence, it could be appropriate to either utilize these two functions or modify storm_PW to take those data into consideration.","category":"page"},{"location":"MoistureMaximization/#Maximization-ratio-and-PMP-estimation","page":"PMP estimation by moisture maximization","title":"Maximization ratio and PMP estimation","text":"","category":"section"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"To maximize the PMP magnitude storms, the precipitation has to be multiplied by the maximization ratio (MR), which is defined by fracPW_maxPW_storm. PW_max can be the maximum precipitable water recorded during the month of the event of the entire dataset or the return value for a given return period for the month of the event. We calculate both fot all the months of interest (May to October): ","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"# pw dataset construction\npw_dataset = DataFrame(Date = dewpoint.DateTime)\npw_dataset.PW = dewpoint_to_PW.(dewpoint.Dew)\n\npw_max = PW_max(pw_dataset.PW, pw_dataset.Date)                  # PW max\npw_100 = PW_return_period(pw_dataset.PW, pw_dataset.Date, 100)   # PW 100\n\nprintln(\"\") # hide","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"If we want information (date, maximization ratio, effective precipitation (EP) and maximized rain) on all the maximized storms, we can call the storm_maximization function : ","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"maximized_storm_PW_max = storm_maximization(storm.Rain, pw_storm.PW, storm.Date, pw_max)\nmaximized_storm_PW_100 = storm_maximization(storm.Rain, pw_storm.PW, storm.Date, pw_100)\nprintln(\"\") # hide","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"Alternatively, we could exclusively estimate the PMP as follows :","category":"page"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"PMP_max = PMP_mm(storm.Rain, pw_storm.PW, storm.Date, pw_max)\nPMP_100 = PMP_mm(storm.Rain, pw_storm.PW, storm.Date, pw_100)\nprintln(\"\") # hide","category":"page"},{"location":"GettingStarted/#Getting-started","page":"Getting started","title":"Getting started","text":"","category":"section"},{"location":"GettingStarted/","page":"Getting started","title":"Getting started","text":"Blah blah blah.","category":"page"},{"location":"#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"Modules =[PMP]","category":"page"},{"location":"#PMP.PearsonType1","page":"Index","title":"PMP.PearsonType1","text":"PearsonType1(a,b,α,β)\n\nThe Pearson Type 1 distribution with shape parameters α and β defined on the interval (a, b) has the probability density function for ayb\n\nf(y a b alpha beta) = frac1B(alpha beta) frac(y-a)^alpha-1 (b-y)^beta-1(b-a)^alpha+beta-1\n\nand 0 elsewhere.\n\nPearsonType1()   # Pearson Type 1 distribution on the unit interval with shape parameters (1,1) i.e. Uniform(0, 1)\n\n\nparams(d)        # Get the parameters, i.e. (a, b, α, β)\n\nlocation(d)      # Get the location parameter, i.e. a\nscale(d)         # Get the scale parameter, i.e. b-a\nshape(d)         # Get the shape parameters, i.e. (α, β)\n\nExternal links\n\nPearson Type 1 distribution on Wikipedia\n\n\n\n\n\n","category":"type"},{"location":"#PMP.PMP_GEV-Tuple{Vector{<:Real}, Vector{Dates.DateTime}, Real, Int64, Int64}","page":"Index","title":"PMP.PMP_GEV","text":"PMP_GEV(rain::Vector{<:Real}, date::Vector{DateTime}, return_time::Real, d₁::Int, d₂::Int)\nPMP_GEV(rain::Vector{<:Real}, date::Vector{DateTime}, return_time::Real)\n\nEstimation of the PMP by univariate GEV method, with a chosen return time.\n\n\n\n\n\n","category":"method"},{"location":"#PMP.PMP_Hershfield-Tuple{Vector{<:Real}, Vector{Dates.DateTime}, Real, Int64, Int64}","page":"Index","title":"PMP.PMP_Hershfield","text":"PMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{DateTime}, K::Real, d₁::Int, d₂::Int)\nPMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{DateTime}, K::Real)\nPMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{DateTime}, d₁::Int, d₂::Int)\nPMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{DateTime})\n\nEstimation of the PMP by Hershfield empirical method.\n\n\n\n\n\n","category":"method"},{"location":"#PMP.PMP_mm","page":"Index","title":"PMP.PMP_mm","text":"PMP_mm(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{DateTime}, pw_max::Vector{<:Real})\nPMP_mm(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{Date}, pw_max::Vector{<:Real})\n\nEstimation of the PMP by moisture maximization.\n\n\n\n\n\n","category":"function"},{"location":"#PMP.PW_max","page":"Index","title":"PMP.PW_max","text":"PW_max(pw::Vector{<:Real}, date::Vector{DateTime})\nPW_max(pw::Vector{<:Real}, date::Vector{Date})\n\nEstimate the maximum precipitable water for each month of interest.\n\n\n\n\n\n","category":"function"},{"location":"#PMP.PW_return_period","page":"Index","title":"PMP.PW_return_period","text":"PW_return_period(pw::Vector{<:Real}, date::Vector{DateTime}, return_period::Int=100)\nPW_return_period(pw::Vector{<:Real}, date::Vector{Date}, return_period::Int=100)\n\nEstimate the precipitable water return value for each month of interest for a given return period.\n\n\n\n\n\n","category":"function"},{"location":"#PMP.dataset-Tuple{String}","page":"Index","title":"PMP.dataset","text":"dataset(name::String)::DataFrame\n\nLoad the dataset associated with name.\n\nDatasets available:\n\nrain: observed precipitations (in mm) recorded at the Montréal-Trudeau International Airport;\ndewpoint: observed dew point (in °C) recorded at the Montréal-Trudeau International Airport.\n\nExamples\n\njulia> PMP.dataset(\"rain\")\n\n\n\n\n\n","category":"method"},{"location":"#PMP.dewpoint_to_PW-Tuple{Real}","page":"Index","title":"PMP.dewpoint_to_PW","text":"dewpoint_to_PW(dew_data::Real)\n\nConvert dew point observation in precipitable water (PW).\n\nThe relation is given by the Table A.1.1 of the annex of the \"Manual on Estimation of Probable Maximum  Precipitation (PMP)\" (WMO, 2009).\n\n\n\n\n\n","category":"method"},{"location":"#PMP.get_max_persisting_dew","page":"Index","title":"PMP.get_max_persisting_dew","text":"get_max_persisting_dew(dewpoint::Vector{<:Real}, frequency::Int, time_int::Int=12)\n\nGet the maximum persisting dewpoint of a storm for which data are taken at a given frequency.\n\nThe highest persisting dewpoint for some specified time interval is the value equalled or exceeded at all observations during the period (WMO, 2009).\n\n\n\n\n\n","category":"function"},{"location":"#PMP.storm_PW","page":"Index","title":"PMP.storm_PW","text":"storm_PW(storm_date::Vector{DateTime}, dewpoint::Vector{<:Real}, dewpoint_date::Vector{DateTime}, d₂::Int, frequency::Int, time_int::Int=12)\nstorm_PW(storm_date::Vector{Date}, dewpoint::Vector{<:Real}, dewpoint_date::Vector{DateTime}, d₂::Int, frequency::Int, time_int::Int=12)\nstorm_PW(storm_date::Vector{DateTime}, dewpoint::Vector{<:Real}, dewpoint_date::Vector{Date}, d₂::Int, frequency::Int, time_int::Int=12)\nstorm_PW(storm_date::Vector{Date}, dewpoint::Vector{<:Real}, dewpoint_date::Vector{Date}, d₂::Int, frequency::Int, time_int::Int=12)\n\nGet the precipitable water for each storm.\n\n\n\n\n\n","category":"function"},{"location":"#PMP.storm_maximization","page":"Index","title":"PMP.storm_maximization","text":"storm_maximization(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{DateTime}, pw_max::Vector{<:Real})\nstorm_maximization(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{Date}, pw_max::Vector{<:Real})\n\nEstimation of the maximization ratio, effective precipitation and maximized precipitation.\n\n\n\n\n\n","category":"function"},{"location":"#PMP.storm_selection","page":"Index","title":"PMP.storm_selection","text":"storm_selection(rain::Vector{<:Real}, date::Vector{DateTime}, p::Real, d₁::Int, d₂::Int=72)\nstorm_selection(rain::Vector{<:Real}, date::Vector{Date}, p::Real, d₁::Int, d₂::Int=72)\nstorm_selection(rain::Vector{<:Real}, date::Vector{DateTime}, p::Real)\nstorm_selection(rain::Vector{<:Real}, date::Vector{Date}, p::Real)\n\nSelect storms of duration d₂ to be maximized from data of duration d₁. \n\n\n\n\n\n","category":"function"},{"location":"#PMP.total_precipitation","page":"Index","title":"PMP.total_precipitation","text":"total_precipitation(rain::Vector{<:Real}, date::Vector{DateTime}, d₁::Int, d₂::Int)\ntotal_precipitation(rain::Vector{<:Real}, date::Vector{Date}, d₁::Int, d₂::Int)\n\nEstimate the greatest precipitation taken over a given duration d₁ on a longer duration d₂.\n\nThe function choose the greatest precipitation of duration d₂ while avoiding taking the same observation twice. Datasets rain and date should not contain missing data. \n\n\n\n\n\n","category":"function"}]
}
