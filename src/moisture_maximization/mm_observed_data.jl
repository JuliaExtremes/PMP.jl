# Storm selection
"""
Selection of storms to be maximized
"""
function storm_selection_cluster(rain_daily::Vector{<:Real}, date::Vector{Dates.Date}, u::Real)

    cluster = getcluster(rain_daily, u, 1)

    storm_rain = DataFrame()
    for i in eachindex(cluster)
        event = DataFrame(Date = date[cluster[i].position[1]], Rain_max = maximum(cluster[i].value), 
                            Rain_total = sum(cluster[i].value), Duration = length(cluster[i].position))
        append!(storm_rain, event)
    end

    return storm_rain
end

function storm_selection_fixed(rain_daily::Vector{<:Real}, date::Vector{Dates.Date}, u::Real)
    
    df = DataFrame(Date = date, Rain = rain_daily)
    df = filter(r -> r.Rain > u, df) 

    # est-ce que le declustering est nécessaire ?
    storm_rain = DataFrame()
    i = 1
    n = size(df, 1)
    while i < n 
        date_dif = Dates.value.(df.Date[i+1:end] .- df.Date[i:end-1])
        ind = findfirst(r -> r > 1, date_dif) 
        if typeof(ind) == Nothing
            ind = n - i + 1
        end
        event = filter(r -> r.Rain == maximum(df.Rain[i:i+ind-1]), df[i:i+ind-1, :])[1, :]
        push!(storm_rain, event)
        i = i + ind
    end
    
    return storm_rain
end



# Maximum persisting dewpoint
"""
The highest persisting dewpoint for some specified time interval (generally 12h or 24h) is the value equalled or exceeded at all 
observations during the period. (2009, WMO)
The function takes hourly dew point of a storm and return the max persisting dew point of it.
"""
function get_max_persisting_dew(time_int::Real, dew_hourly::Vector{<:Real})

    persisting_dews = []
    for k = 1:length(dew_hourly) - (time_int-1)
        push!(persisting_dews, minimum(dew_hourly[k:k + (time_int-1)]))
    end

    return maximum(persisting_dews)
end # fonction pour une tempête (dew_hourly sont les données d'une tempête)



# Conversion dewpoint to PW
"""
Conversion from dew point to precipitable water (PW), using the relation given by the Table A.1.1 of the 
annexe of the "Manual on Estimation of Probable Maximum Precipitation (PMP)" (2009, WMO). The function takes 
a single dew point observation as an argument and returns its associated precipitable water.
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
Estimation of the monthly maximum precipitable water 
"""
function PW_max(pw_storm::Vector{<:Real}, date::Vector{Dates.Date}) # nécéssaire ?
    
    month = Dates.month.(date)
    df = DataFrame(PW = pw_storm, Month = month)
    PW_max = combine(groupby(df, :Month), :PW => maximum => :PW_max)
    
    return PW_max
end


"""
Estimation of the monthly precipitable water for given a return period (?)
"""
function PW_return_period(return_period::Real, pw_storm::Vector{<:Real}, date::Vector{Dates.Date})
    
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



# PMP estimation
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