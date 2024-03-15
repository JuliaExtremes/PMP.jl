"""
    PearsonType1c(b,μ,ν)

The *Pearson Type 1 distribution* with shape parameters `μ` and `ν` defined on the interval (0, `b`) has the probability density function for ``0<y<b``
```math
f(y; b, \\mu, \\nu) = \\frac{1}{B(\\mu, \\nu)} \\frac{y^{\\mu\\nu-1} (b-y)^{\\nu(1-\\mu) - 1}}{b^{\\nu-1}},
```
and 0 elsewhere.


```julia
PearsonType1c()   # Pearson Type 1 distribution on the unit interval with shape parameters (1/2,2) i.e. Uniform(0, 1)


params(d)        # Get the parameters, i.e. (b, μ, ν)

scale(d)         # Get the scale parameter, i.e. b
shape(d)         # Get the shape parameters, i.e. (μ, ν)
```

"""
struct PearsonType1c{T<:Real} <: ContinuousUnivariateDistribution
    b::T
    μ::T
    ν::T
    PearsonType1c{T}(b::T, μ::T, ν::T) where {T<:Real} = new{T}(b, μ, ν)
end

function PearsonType1c(b::T, μ::T, ν::T ; check_args::Bool=true) where {T <: Real}
    @check_args PearsonType1c (b, b>0) (μ, μ>zero(μ)) (ν, ν>zero(ν))
    return PearsonType1c{T}(b, μ, ν)
end

PearsonType1c(b::Real, μ::Real, ν::Real; check_args::Bool=true) = PearsonType1c(promote(b, μ, ν)...; check_args=check_args)
PearsonType1c(b::Real, μ::Integer, ν::Integer; check_args::Bool=true) = PearsonType1c(float(b), float(μ), float(ν); check_args=check_args)
PearsonType1c() = PearsonType1c{Float64}(1.0, 0.5, 2.0)



# Parameters
"""
    location(pd::PearsonType1c)

Obtain distribution location parameter 0.
"""
function location(pd::PearsonType1c)
    return 0    
end

"""
    scale(pd::PearsonType1c)

Obtain distribution scale parameter b.
"""
function scale(pd::PearsonType1c)
    return pd.b
end

"""
    shape(pd::PearsonType1c)

Obtain distribution shape parameters μ and ν.
"""
function shape(pd::PearsonType1c)
    return (pd.μ, pd.ν)
end

"""
    shape(pd::PearsonType1c)

Obtain all distribution parameters b, μ and ν.
"""
function params(pd::PearsonType1c)
 return pd.b, shape(pd)...
end



"""
    logpdf(pd::PearsonType1c, x::Real)

Compute the log of the value of the probability density function of pd at point x.
"""
function logpdf(pd::PearsonType1c, x::Real)
    B = loggamma(pd.ν) - (loggamma(pd.μ*pd.ν) + loggamma(pd.ν*(1 - pd.μ)))
    return B + (pd.μ*pd.ν - 1)*log(x) + (pd.ν*(1 - pd.μ) - 1)*log(pd.b - x) - (pd.ν - 1)*log(pd.b)
end



# MLE
"""
    fit_mle(pd::Type{<:PearsonType1c}, y::Vector{<:Real}, initialvalues::Vector{<:Real})
    fit_mle(pd::Type{<:PearsonType1c}, y::Vector{<:Real})

Estimate parameters of a PearsonType1c distribution with likelihood maximization.
"""

function fit_mle(pd::Type{<:PearsonType1c}, y::Vector{<:Real}, initialvalues::Vector{<:Real})
 
    iv = [initialvalues[1], initialvalues[2], initialvalues[3]]
    initialvalues[1] = max(initialvalues[1], maximum(y)) + .01*maximum(y)

    loglike(θ::Vector{<:Real}) = sum(logpdf.(PearsonType1c(θ...),y))
    fobj(θ) = -loglike(θ)

    lower = [maximum(y), 2*eps(), 2*eps()]
    upper = [Inf, Inf, Inf]
    
    res = optimize(fobj, lower, upper, initialvalues, autodiff = :forward)
    
    if Optim.converged(res)
        θ̂ = Optim.minimizer(res)
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        θ̂ = iv
    end
        
    return PearsonType1c(θ̂...)
end

function fit_mle(pd::Type{<:PearsonType1c}, y::Vector{<:Real})
    initialvalues =  getinitialvalues(PearsonType1c, y)
    fd = fit_mle(pd, y, initialvalues)
    return fd
end 



# Mixed method
"""
    getinitialvalues(pd::Type{<:PearsonType1c}, y::Vector{<:Real})

Find initial values for the likelihood maximization algorithm. 

Initial values of shape parameters μ and ν are estimated by method of moments and initial value of scale parameter b is estimated by likelihood maximization.
"""

function getinitialvalues(pd::Type{<:PearsonType1c}, y::Vector{<:Real})
    initialvalues = getinitialvalues(PearsonType1b, y)
    
    α = initialvalues[1]
    β = initialvalues[2]
    b = initialvalues[3]
    
    ν = α + β
    μ = α/ν

    return [b, μ, ν]
end

