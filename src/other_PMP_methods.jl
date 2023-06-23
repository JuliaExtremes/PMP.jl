"""
    PMP_GEV(rain_daily::Vector{<:Real}, date::Vector{Dates.Date}, return_time::Real)

Estimation of the PMP by univariate GEV method, with a chosen return time.
"""
function PMP_GEV(rain_daily::Vector{<:Real}, date::Vector{Dates.Date}, return_time::Real)
    
    df = DataFrame(Rain = rain_daily, Year = Dates.year.(date))
    AMP = combine(groupby(df, :Year), :Rain => maximum => :AMP)

    return returnlevel(gevfit(AMP, :AMP), return_time).value[]

end



"""
    PMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{Dates.Date}, K::Real)

Estimation of the PMP by Hershfield empirical method.
"""
function PMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{Dates.Date}, K::Real)
    
    df = DataFrame(Rain = rain_daily, Year = year.(date))
    AMP = combine(groupby(df, :Year), :Rain => maximum => :AMP)

    return mean(AMP.AMP) + K*std(AMP.AMP)

end 

function PMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{Dates.Date})

    df = DataFrame(Rain = rain_daily, Year = year.(date))
    AMP = combine(groupby(df, :Year), :Rain => maximum => :AMP)
    
    indmax = findmax(AMP.AMP)[2]
    AMP2 = AMP.AMP[Not(indmax)]
    K = (maximum(AMP.AMP) - mean(AMP2))/std(AMP2)

    return mean(AMP.AMP) + K*std(AMP.AMP), K

end 