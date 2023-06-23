"""
    total_precipitation(rain::Vector{<:Real}, date::Vector{Dates.DateTime}, d₁::Int, d₂::Int)

Estimate the greatest precipitation taken over a given duration `d₁` on a longer duration `d₂`.

`rain` and `date` should not contain missing data.
"""
function total_precipitation(rain::Vector{<:Real}, date::Vector{Dates.DateTime}, d₁::Int, d₂::Int)

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
    n = length(df1.Rain)
    
    while i < n
        
        df1.date_dif = Hour.(abs.(storm.Date[i] .- df1.Date))
        
        filter!(r -> r.date_dif >= Hour(d₂), df1)
        select!(df1, :Date, :Rain)
        
        storm = storm[1:i, :]
        append!(storm, df1)
        
        n = length(storm.Rain)
        i += 1
        
        df1 = storm[i:end, :]

    end   
    sort!(storm)
    
    return storm

end

total_precipitation(rain::Vector{<:Real}, date::Vector{Dates.Date}, d₁::Int, d₂::Int) = total_precipitation(rain, Dates.DateTime.(date), d₁::Int, d₂::Int)



# Storm selection
"""
    storm_selection(rain::Vector{<:Real}, date::Vector{Dates.DateTime}, p::Real, d₁::Int, d₂::Int=72)

Select storms to be maximized. 
"""
function storm_selection(rain::Vector{<:Real}, date::Vector{Dates.DateTime}, p::Real, d₁::Int, d₂::Int=72)

    @assert 0 < p <= 1 "p should lie in the interval (0, 1]."
    
    df1 = PMP.total_precipitation(rain, date, d₁, d₂)
    df1.Year = Dates.year.(df1.Date)
    
    nYear = df1.Year[end] - df1.Year[1] + 1
    storm = DataFrame() 
    
    for i in 1:nYear
        storm_year = groupby(df1, :Year)[i]
        nStormMax = floor(Int, p*length(storm_year.Rain)) + 1
        
        storm_max = sort(storm_year, :Rain, rev=true)[1:nStormMax, :]
        append!(storm, storm_max)
    end
    select!(storm, :Date, :Rain)
    storm.Date = Date.(storm.Date)
    
    return storm

end 

storm_selection(rain::Vector{<:Real}, date::Vector{Dates.Date}, p::Real, d₁::Int, d₂::Int=72) = storm_selection(rain, Dates.DateTime.(date), p, d₁, d₂)



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
    PW_max(pw_storm::Vector{<:Real}, date::Vector{Dates.Date})

Estimate the maximum precipitable water for each month of interest.
"""
function PW_max(pw_storm::Vector{<:Real}, date::Vector{Dates.Date})
    
    month = Dates.month.(date)
    df = DataFrame(PW = pw_storm, Month = month)
    PW_max = combine(groupby(df, :Month), :PW => maximum => :PW_max)
    
    return PW_max
end


"""
    PW_return_period(pw_storm::Vector{<:Real}, date::Vector{Dates.Date}, return_period::Int=100)

Estimate the precipitable water return value for each month of interest.
"""
function PW_return_period(pw_storm::Vector{<:Real}, date::Vector{Dates.Date}, return_period::Int=100)
    
    ym = Dates.yearmonth.(date)
    month = Dates.month.(date)
    df = DataFrame(PW = pw_storm, YM = ym, Month = month)
    PW_month = combine(groupby(df, :YM), :PW => maximum => :PW, :Month => :Month)
    ind = combine(groupby(df, :Month), :PW => maximum => :PW)
    
    PW_rp = DataFrame()
    for i in ind.Month
        df = filter(:Month => ==(i), PW_month)
        pw_rp_month = DataFrame(Month = i, PW_rp = returnlevel(gevfit(df, :PW), return_period).value[1])
        append!(PW_rp, pw_rp_month)
    end 
# x1.2 
    return PW_rp
end



# Storm maximization
"""
Estimation of the maximization ratio, effective precipitation, maximized storm and PMP (moisture maximization method).
"""
function PMP_mm(rain_storm::Vector{<:Real}, pw_storm::Vector{<:Real}, date_storm::Vector{Dates.Date}, pw_max::Vector{<:Real})
    months = Dates.month.(date_storm) 
    storm = DataFrame(Rain = rain_storm, PW = pw_storm, PW_max = pw_max[months .- (minimum(months)-1)])

    storm.EP = storm.Rain ./ storm.PW

    storm.maximization_ratio = storm.PW_max ./ storm.PW
    storm.bounded_maximization_ratio = min.(2, storm.maximization_ratio)

    storm.maximized_rain = storm.Rain .* storm.maximization_ratio
    storm.bounded_maximized_rain = storm.Rain .* storm.bounded_maximization_ratio

    return maximum(storm.maximized_rain), maximum(storm.bounded_maximized_rain), storm
end

# j'aimerais diviser cette fonction en deux : une fonction de maximisation des tempetes qui retourne un dataframe 
# contenant toutes les tempetes maximisees et une seconde qui retourne la PMP (devrait-on seulement retourner 
# l'evenement associe a la PMP sans toutes les tempetes ?)
# modifier pour le ratio de maximisation borne 