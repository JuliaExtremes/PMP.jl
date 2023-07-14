"""
    PMP_GEV(rain::Vector{<:Real}, date::Vector{DateTime}, return_time::Real, d₁::Int, d₂::Int)
    PMP_GEV(rain::Vector{<:Real}, date::Vector{DateTime}, return_time::Real)

Estimation of the PMP by univariate GEV method, with a chosen return time.
"""
function PMP_GEV(rain::Vector{<:Real}, date::Vector{DateTime}, return_time::Real, d₁::Int, d₂::Int)
    
    df = PMP.total_precipitation(rain, date, d₁, d₂)
    df.Year = year.(df.Date)
    AMP = combine(groupby(df, :Year), :Rain => maximum => :AMP)

    return returnlevel(gevfit(AMP, :AMP), return_time).value[]

end

function PMP_GEV(rain::Vector{<:Real}, date::Vector{DateTime}, return_time::Real)
    
    df = DataFrame(Rain = rain, Year = year.(date))
    AMP = combine(groupby(df, :Year), :Rain => maximum => :AMP)

    return returnlevel(gevfit(AMP, :AMP), return_time).value[]

end

PMP_GEV(rain::Vector{<:Real}, date::Vector{Date}, return_time::Real, d₁::Int, d₂::Int) = PMP_GEV(rain, DateTime.(date), return_time, d₁, d₂)
PMP_GEV(rain::Vector{<:Real}, date::Vector{Date}, return_time::Real) = PMP_GEV(rain, DateTime.(date), return_time)



# ajuster selon la taille de l'échantillon p.66
"""
    PMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{DateTime}, K::Real, d₁::Int, d₂::Int)
    PMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{DateTime}, K::Real)
    PMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{DateTime}, d₁::Int, d₂::Int)
    PMP_Hershfield(rain_daily::Vector{<:Real}, date::Vector{DateTime})

Estimation of the PMP by Hershfield empirical method.
"""
function PMP_Hershfield(rain::Vector{<:Real}, date::Vector{DateTime}, K::Real, d₁::Int, d₂::Int)
    
    df = PMP.total_precipitation(rain, date, d₁, d₂)
    df.Year = year.(df.Date)
    AMP = combine(groupby(df, :Year), :Rain => maximum => :AMP)

    return mean(AMP.AMP) + K*std(AMP.AMP)

end 

function PMP_Hershfield(rain::Vector{<:Real}, date::Vector{DateTime}, K::Real)
    
    df = DataFrame(Rain = rain, Year = year.(date))
    AMP = combine(groupby(df, :Year), :Rain => maximum => :AMP)

    return mean(AMP.AMP) + K*std(AMP.AMP)

end 

function PMP_Hershfield(rain::Vector{<:Real}, date::Vector{DateTime}, d₁::Int, d₂::Int)

    df = PMP.total_precipitation(rain, date, d₁, d₂)
    df.Year = year.(df.Date)
    AMP = combine(groupby(df, :Year), :Rain => maximum => :AMP)
    
    indmax = findmax(AMP.AMP)[2]
    AMP2 = AMP.AMP[Not(indmax)]
    K = (maximum(AMP.AMP) - mean(AMP2))/std(AMP2)

    return mean(AMP.AMP) + K*std(AMP.AMP), K

end 

function PMP_Hershfield(rain::Vector{<:Real}, date::Vector{DateTime})

    df = DataFrame(Rain = rain, Year = year.(date))
    AMP = combine(groupby(df, :Year), :Rain => maximum => :AMP)
    
    indmax = findmax(AMP.AMP)[2]
    AMP2 = AMP.AMP[Not(indmax)]
    K = (maximum(AMP.AMP) - mean(AMP2))/std(AMP2)

    return mean(AMP.AMP) + K*std(AMP.AMP), K

end 

PMP_Hershfield(rain::Vector{<:Real}, date::Vector{Date}, K::Real, d₁::Int, d₂::Int) = PMP_Hershfield(rain, DateTime.(date), K, d₁, d₂)
PMP_Hershfield(rain::Vector{<:Real}, date::Vector{Date}, K::Real) = PMP_Hershfield(rain, DateTime.(date), K)
PMP_Hershfield(rain::Vector{<:Real}, date::Vector{Date}, d₁::Int, d₂::Int) = PMP_Hershfield(rain, DateTime.(date), d₁, d₂)
PMP_Hershfield(rain::Vector{<:Real}, date::Vector{Date}) = PMP_Hershfield(rain, DateTime.(date))