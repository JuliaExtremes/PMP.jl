"""
    total_precipitation(rain::Vector{<:Real}, date::Vector{DateTime}, d₁::Int, d₂::Int)

Estimate the greatest precipitation taken over a given duration `d₁` on a longer duration `d₂`.

`rain` and `date` should not contain missing data.
"""
function total_precipitation(rain::Vector{<:Real}, date::Vector{DateTime}, d₁::Int, d₂::Int)

    @assert d₁ <= d₂ "the second duration should be longer than the first one."
    @assert length(rain) == length(date) "the vectors of rain data and the associated dates should be of the same length."

    nb = floor(Int, d₂/d₁) # window size

    rain_d₂ = RollingFunctions.rolling(sum, rain, nb)
    size = length(rain_d₂)

    df1 = DataFrame(Rain = rain_d₂, Date = date[1:size])
    filter!(r -> r.Rain>0, df1)
    sort!(df1, rev=true)
    df1 = combine(groupby(df1, :Date), :Rain => maximum => :Rain)

    storm = DataFrame(df1[1, :])
    
    i = 1
    n = nrow(df1)
    
    while i < n
        
        df1.date_dif = Hour.(abs.(storm.Date[i] .- df1.Date))
        
        filter!(r -> r.date_dif >= Hour(d₂), df1)
        select!(df1, :Date, :Rain)
        
        storm = storm[1:i, :]
        append!(storm, df1)
        
        n = nrow(storm)
        i += 1
        
        df1 = storm[i:end, :]

    end   
    sort!(storm)
    
    return storm

end

total_precipitation(rain::Vector{<:Real}, date::Vector{Date}, d₁::Int, d₂::Int) = total_precipitation(rain, DateTime.(date), d₁::Int, d₂::Int)



# Storm selection
"""
    storm_selection(rain::Vector{<:Real}, date::Vector{DateTime}, p::Real, d₁::Int, d₂::Int=72)

Select storms to be maximized. 
"""
function storm_selection(rain::Vector{<:Real}, date::Vector{DateTime}, p::Real, d₁::Int, d₂::Int=72)

    @assert 0 < p <= 1 "p should lie in the interval (0, 1]."
    
    df1 = PMP.total_precipitation(rain, date, d₁, d₂)
    df1.Year = year.(df1.Date)
    
    nYear = df1.Year[end] - df1.Year[1] + 1
    storm_year = groupby(df1, :Year)
    storm = DataFrame() 
    
    for i in 1:nYear

        nStormMax = floor(Int, p*nrow(storm_year[i])) + 1
        storm_max = sort(storm_year[i], :Rain, rev=true)[1:nStormMax, :]
        append!(storm, storm_max)

    end
    
    select!(storm, :Date, :Rain)
    storm.Date = Date.(storm.Date)
    
    return storm

end 

function storm_selection(rain::Vector{<:Real}, date::Vector{DateTime}, p::Real)

    @assert 0 < p <= 1 "p should lie in the interval (0, 1]."
    
    df1 = DataFrame(Rain = rain, Date = date, Year = year.(date))
    
    nYear = df1.Year[end] - df1.Year[1] + 1
    storm_year = groupby(df1, :Year)
    storm = DataFrame() 
    
    for i in 1:nYear

        nStormMax = floor(Int, p*nrow(storm_year[i])) + 1
        storm_max = sort(storm_year[i], :Rain, rev=true)[1:nStormMax, :]
        append!(storm, storm_max)

    end
    
    select!(storm, :Date, :Rain)
    storm.Date = Date.(storm.Date)
    
    return storm

end

storm_selection(rain::Vector{<:Real}, date::Vector{Date}, p::Real, d₁::Int, d₂::Int=72) = storm_selection(rain, DateTime.(date), p, d₁, d₂)
storm_selection(rain::Vector{<:Real}, date::Vector{Date}, p::Real) = storm_selection(rain, DateTime.(date), p)



# Maximum persisting dewpoint
"""
    get_max_persisting_dew(dew_hourly::Vector{<:Real}, frequency::Int, time_int::Int=12)

Get the maximum persisting dew point of a storm for which data are taken at a given frequency.
    
The highest persisting dewpoint for some specified time interval is the value equalled or exceeded at all 
observations during the period (2009, WMO).
"""
function get_max_persisting_dew(dew_hourly::Vector{<:Real}, frequency::Int, time_int::Int=12)

    nb = floor(Int, time_int/frequency) # window size

    persisting_dews = maximum(RollingFunctions.rollmin(dew_hourly, nb))

    return persisting_dews

end



# Dew point to PW
"""
    dewpoint_to_PW(dew_data::Real)

Convert dew point observation in precipitable water (PW).

The relation is given by the Table A.1.1 of the annex of the "Manual on Estimation of Probable Maximum 
Precipitation (PMP)" (2009, WMO).
"""
function dewpoint_to_PW(dew_data::Real)

    pw = [8., 9, 10, 11, 12, 13, 15, 16, 18, 19, 21, 23, 25, 28, 30, 33, 36, 40, 44, 48, 
    52, 57, 62, 68, 74, 81, 88, 96, 105, 114, 123] 
    nonogram = DataFrame(dewpoint = collect(0.:1:30), pw = pw)

    if dew_data < 0.
        return 8.
    end
    if dew_data > 30.
        return 123.
    end

    lower_point = filter(r-> r.dewpoint == floor(dew_data), nonogram).pw[1]
    upper_point = filter(r-> r.dewpoint == floor(dew_data) + 1., nonogram).pw[1]
    
    pw_data = (dew_data - floor(dew_data)) * (upper_point - lower_point) + lower_point
    
    return pw_data

end



# Maximum precipitable water 
"""
    PW_max(pw_storm::Vector{<:Real}, date::Vector{DateTime})

Estimate the maximum precipitable water for each month of interest.
"""
function PW_max(pw_storm::Vector{<:Real}, date::Vector{DateTime})
    
    month = Dates.month.(date)
    df = DataFrame(PW = pw_storm, Month = month)
    PW_max = combine(groupby(df, :Month), :PW => maximum => :PW_max)
    
    return PW_max

end

PW_max(pw_storm::Vector{<:Real}, date::Vector{Date}) = PW_max(pw_storm, DateTime.(date))



"""
    PW_return_period(pw_storm::Vector{<:Real}, date::Vector{DateTime}, return_period::Int=100)

Estimate the precipitable water return value for each month of interest.
"""
function PW_return_period(pw_storm::Vector{<:Real}, date::Vector{DateTime}, return_period::Int=100)
    
    ym = Dates.yearmonth.(date)
    month = Dates.month.(date)
    df1 = DataFrame(PW = pw_storm, YM = ym, Month = month)
    
    PW_month = combine(groupby(df1, :YM), :PW => maximum => :PW, :Month => :Month)
    ind = combine(groupby(df1, :Month), :PW => maximum => :PW)
    
    PW_rp = DataFrame(Month = Int[], PW_rp = Float64[])
    
    for i in ind.Month

        df2 = filter(:Month => ==(i), PW_month)
        pw_rp_month = DataFrame(Month = i, PW_rp = returnlevel(gevfit(df2, :PW), return_period).value)
        append!(PW_rp, pw_rp_month)

    end 
 
    return PW_rp
end

PW_return_period(pw_storm::Vector{<:Real}, date::Vector{Date}, return_period::Int=100) = PW_return_period(pw_storm, DateTime.(date), return_period)



# Storm maximization
"""
    storm_maximization(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{DateTime}, pw_max::Vector{<:Real})

Estimation of the maximization ratio, effective precipitation and maximized precipitation.
"""
function storm_maximization(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{DateTime}, pw_max::Vector{<:Real})

    months = Dates.month.(date_storm)
    storm = DataFrame(Rain = rain_storm, PW = pw_storm, PW_max = pw_max[months .- (minimum(months)-1)])

    EP = rain_storm ./ pw_storm
    maximization_ratio = pw_max[months .- (minimum(months)-1)] ./ pw_storm

    maximized_rain = rain_storm .* maximization_ratio

    storm = DataFrame(Date = date_storm, MR = maximization_ratio, EP = EP, Maximized_Rain = maximized_rain)

    return storm

end

storm_maximization(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{Date}, pw_max::Vector{<:Real}) = storm_maximization(rain_storm, pw_storm, DateTime.(date_storm), pw_max)



# PMP_mm
"""
    PMP_mm(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{DateTime}, pw_max::Vector{<:Real})

Estimation of the PMP by moisture maximization.
"""
function PMP_mm(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{DateTime}, pw_max::Vector{<:Real})
    
    storm = PMP.storm_maximization(rain_storm, pw_storm, date_storm, pw_max)

    return maximum(storm.Maximized_Rain)

end

PMP_mm(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{Date}, pw_max::Vector{<:Real}) = PMP_mm(rain_storm, pw_storm, DateTime.(date_storm), pw_max)