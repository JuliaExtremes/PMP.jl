var documenterSearchIndex = {"docs":
[{"location":"MoistureMaximization/#PMP-estimation-by-moisture-maximization","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"","category":"section"},{"location":"MoistureMaximization/","page":"PMP estimation by moisture maximization","title":"PMP estimation by moisture maximization","text":"Blah blah blah","category":"page"},{"location":"MoistureMaximization/#Subsection-1","page":"PMP estimation by moisture maximization","title":"Subsection 1","text":"","category":"section"},{"location":"MoistureMaximization/#Subsection-2","page":"PMP estimation by moisture maximization","title":"Subsection 2","text":"","category":"section"},{"location":"GettingStarted/#Getting-started","page":"Getting started","title":"Getting started","text":"","category":"section"},{"location":"GettingStarted/","page":"Getting started","title":"Getting started","text":"Blah blah blah.","category":"page"},{"location":"#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"Modules =[PMP]","category":"page"},{"location":"#PMP.PearsonType1","page":"Index","title":"PMP.PearsonType1","text":"PearsonType1(a,b,α,β)\n\nThe Pearson Type 1 distribution with shape parameters α and β defined on the interval (a, b) has the probability density function for ayb\n\nf(y a b alpha beta) = frac1B(alpha beta) frac(y-a)^alpha-1 (b-y)^beta-1(b-a)^alpha+beta-1\n\nand 0 elsewhere.\n\nPearsonType1()   # Pearson Type 1 distribution on the unit interval with shape parameters (1,1) i.e. Uniform(0, 1)\n\n\nparams(d)        # Get the parameters, i.e. (a, b, α, β)\n\nlocation(d)      # Get the location parameter, i.e. a\nscale(d)         # Get the scale parameter, i.e. b-a\nshape(d)         # Get the shape parameters, i.e. (α, β)\n\nExternal links\n\nPearson Type 1 distribution on Wikipedia\n\n\n\n\n\n","category":"type"},{"location":"#PMP.PMP_GEV-Tuple{Vector{<:Real}, Vector{Dates.DateTime}, Real, Int64, Int64}","page":"Index","title":"PMP.PMP_GEV","text":"PMP_GEV(rain::Vector{<:Real}, date::Vector{DateTime}, return_time::Real, d₁::Int, d₂::Int)\nPMP_GEV(rain::Vector{<:Real}, date::Vector{DateTime}, return_time::Real)\n\nEstimation of the PMP by univariate GEV method, with a chosen return time.\n\n\n\n\n\n","category":"method"},{"location":"#PMP.PMP_Hershfield-Tuple{Vector{<:Real}, Vector{Dates.DateTime}, Real, Int64, Int64}","page":"Index","title":"PMP.PMP_Hershfield","text":"PMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{DateTime}, K::Real, d₁::Int, d₂::Int)\nPMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{DateTime}, K::Real)\nPMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{DateTime}, d₁::Int, d₂::Int)\nPMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{DateTime})\n\nEstimation of the PMP by Hershfield empirical method.\n\n\n\n\n\n","category":"method"},{"location":"#PMP.PMP_mm","page":"Index","title":"PMP.PMP_mm","text":"PMP_mm(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{DateTime}, pw_max::Vector{<:Real})\nPMP_mm(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{Date}, pw_max::Vector{<:Real})\n\nEstimation of the PMP by moisture maximization.\n\n\n\n\n\n","category":"function"},{"location":"#PMP.PW_max","page":"Index","title":"PMP.PW_max","text":"PW_max(pw_storm::Vector{<:Real}, date::Vector{DateTime})\nPW_max(pw_storm::Vector{<:Real}, date::Vector{Date})\n\nEstimate the maximum precipitable water for each month of interest.\n\n\n\n\n\n","category":"function"},{"location":"#PMP.PW_return_period","page":"Index","title":"PMP.PW_return_period","text":"PW_return_period(pw_storm::Vector{<:Real}, date::Vector{DateTime}, return_period::Int=100)\nPW_return_period(pw_storm::Vector{<:Real}, date::Vector{Date}, return_period::Int=100)\n\nEstimate the precipitable water return value for each month of interest for a given return period.\n\n\n\n\n\n","category":"function"},{"location":"#PMP.dewpoint_to_PW-Tuple{Real}","page":"Index","title":"PMP.dewpoint_to_PW","text":"dewpoint_to_PW(dew_data::Real)\n\nConvert dew point observation in precipitable water (PW).\n\nThe relation is given by the Table A.1.1 of the annex of the \"Manual on Estimation of Probable Maximum  Precipitation (PMP)\" (WMO, 2009).\n\n\n\n\n\n","category":"method"},{"location":"#PMP.get_max_persisting_dew","page":"Index","title":"PMP.get_max_persisting_dew","text":"get_max_persisting_dew(dew_hourly::Vector{<:Real}, frequency::Int, time_int::Int=12)\n\nGet the maximum persisting dew point of a storm for which data are taken at a given frequency.\n\nThe highest persisting dewpoint for some specified time interval is the value equalled or exceeded at all observations during the period (WMO, 2009).\n\n\n\n\n\n","category":"function"},{"location":"#PMP.storm_maximization","page":"Index","title":"PMP.storm_maximization","text":"storm_maximization(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{DateTime}, pw_max::Vector{<:Real})\nstorm_maximization(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{Date}, pw_max::Vector{<:Real})\n\nEstimation of the maximization ratio, effective precipitation and maximized precipitation.\n\n\n\n\n\n","category":"function"},{"location":"#PMP.storm_selection","page":"Index","title":"PMP.storm_selection","text":"storm_selection(rain::Vector{<:Real}, date::Vector{DateTime}, p::Real, d₁::Int, d₂::Int=72)\nstorm_selection(rain::Vector{<:Real}, date::Vector{Date}, p::Real, d₁::Int, d₂::Int=72)\nstorm_selection(rain::Vector{<:Real}, date::Vector{DateTime}, p::Real)\nstorm_selection(rain::Vector{<:Real}, date::Vector{Date}, p::Real)\n\nSelect storms of duration d₂ to be maximized from data of duration d₁. \n\n\n\n\n\n","category":"function"},{"location":"#PMP.total_precipitation","page":"Index","title":"PMP.total_precipitation","text":"total_precipitation(rain::Vector{<:Real}, date::Vector{DateTime}, d₁::Int, d₂::Int)\ntotal_precipitation(rain::Vector{<:Real}, date::Vector{Date}, d₁::Int, d₂::Int)\n\nEstimate the greatest precipitation taken over a given duration d₁ on a longer duration d₂.\n\nThe function choose the greatest precipitation of duration d₂ while avoiding taking the same observation twice. Datasets rain and date should not contain missing data. \n\n\n\n\n\n","category":"function"}]
}
