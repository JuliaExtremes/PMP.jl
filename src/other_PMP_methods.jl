"""PMP estimation with GEV method. Return the PMP for the chosen return time, the plot of the maximal annual 
precipitation and the GEV diagnostic plot."""

function PMP_GEV(rain_daily::Vector{<:Real}, date::Vector{Dates.Date}, return_time::Real)
    data = DataFrame(Rain = rain_daily, Year = Dates.year.(date))
    AMP = combine(groupby(data, :Year), :Rain => maximum => :AMP)
    gevfm = gevfit(AMP, :AMP)
    PMP = returnlevel(gevfm, return_time).value[]

    p1 = plot(AMP, x=:Year, y=:AMP, Geom.line, Geom.point, Guide.Title("Maximal Annual Precipitation (mm)"), 
                yintercept=[minimum(AMP.AMP), maximum(AMP.AMP)], Geom.hline(style=:dot))
    p2 = diagnosticplots(gevfm)
    
    return PMP, p1, p2
end



"""Hershfield empirical method"""
function PMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{Dates.Date}, K::Real)
    data = DataFrame(Rain = rain_daily, Year = year.(date))
    AMP = combine(groupby(data, :Year), :Rain => maximum => :AMP)
    
    X̄ = mean(AMP.AMP)
    S = std(AMP.AMP)
    
    PMP = X̄ + K*S

    return PMP
end 

function PMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{Dates.Date})
    data = DataFrame(Rain = rain_daily, Year = year.(date))
    AMP = combine(groupby(data, :Year), :Rain => maximum => :AMP)
    
    X̄ = mean(AMP.AMP)
    S = std(AMP.AMP)
    
    indmax = findmax(AMP.AMP)[2]
    AMP2 = AMP.AMP[Not(indmax)]
    K = (maximum(AMP.AMP) - mean(AMP2))/std(AMP2)
    
    PMP = X̄ + K*S

    return PMP, K
end 