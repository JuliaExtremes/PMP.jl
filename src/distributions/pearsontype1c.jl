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
    @check_args PearsonType1c (b, b>0) (μ, 1>=μ>zero(μ)) (ν, ν>zero(ν))
    return PearsonType1c{T}(b, μ, ν)
end

PearsonType1c(b::Real, μ::Real, ν::Real; check_args::Bool=true) = PearsonType1c(promote(b, μ, ν)...; check_args=check_args)
PearsonType1c(b::Integer, μ::Real, ν::Integer; check_args::Bool=true) = PearsonType1c(float(b), float(μ), float(ν); check_args=check_args)
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



function betalogpdf_reparam(μ::Real, ν::Real, x::Real) # reparam of StatsFuns betalogpdf
    y = clamp(x, 0, 1)
    val = xlogy(μ*ν - 1., y) + xlog1py(ν*(1. - μ) - 1., -y) - logbeta(μ*ν, ν*(1. - μ))
    return x < 0 || x > 1 ? oftype(val, -Inf) : val
end



"""
    logpdf(pd::PearsonType1c, x::Real)

Compute the log of the value of the probability density function of pd at point x.
"""
function logpdf(pd::PearsonType1c, x::Real)
    return betalogpdf_reparam(pd.μ, pd.ν, x/pd.b) - log(pd.b)
end

"""
function logpdf(pd::PearsonType1b, x::Real)
    return logpdf(Beta(shape(pd)...), x/pd.b) - log(pd.b)
end
"""

#function logpdf(pd::PearsonType1c, x::Real)
#    B, sign = logabsbeta(pd.μ * pd.ν, pd.ν * (1. - pd.μ))
#    if sign < 0 
#        throw(DomainError((pd.μ, pd.ν), "`beta(a, b)` must be non-negative"))
#    elseif pd.b - x <= 0
#        throw(DomainError((pd.b, x), "Argument to logpdf must be within the domain of the distribution"))
#    else
#        return -B + (pd.μ * pd.ν - 1.) * log(x) + (pd.ν * (1. - pd.μ) - 1.) * log(pd.b - x) - (pd.ν - 1.) * log(pd.b)
#    end
#end



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
    upper = [Inf, 1, Inf]
    
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
    
    b = initialvalues[1]
    α = initialvalues[2]
    β = initialvalues[3]
    
    ν = α + β
    μ = α/ν

    return [b, μ, ν]
end



# Bayesian fitting
"""
    fit_bayes(pd::Type{<:PearsonType1c}, prior::ContinuousUnivariateDistribution, y::Vector{<:Real}, niter::Int, warmup::Int)

Estimate parameters of a PearsonType1c distribution with Bayesian inference.

Use NUTS (No U-Turn) sampler. The prior refers to the prior distribution of the bound parameter b.
"""

function fit_bayes(pd::Type{<:PearsonType1c}, prior::ContinuousUnivariateDistribution, y::Vector{<:Real}, niter::Int, warmup::Int)
    nparam = 3
    iv = getinitialvalues(pd, y)
    initialvalues = [log(iv[1]), logit(iv[2]), log(iv[3])]

    # Defines the llh function and the gradient for the NUTS algo
    logf(θ::DenseVector) = sum(logpdf.(pd(exp(θ[1]), logistic(θ[2]), exp(θ[3])), y)) + logpdf(prior, exp(θ[1])) 
    Δlogf(θ::DenseVector) = ForwardDiff.gradient(logf, θ)
    function logfgrad(θ::DenseVector)
        ll = logf(θ)
        g = Δlogf(θ)
        return ll, g
    end

    # NUTS algo
    sim = Chains(niter, nparam, start=(warmup+1))
    θ = NUTSVariate(initialvalues, logfgrad)
    for i in 1:niter
        MambaLite.sample!(θ, adapt=(i<=warmup))
        if i>warmup
            sim[i, :, 1] = θ
        end
    end
    
    b̂ = exp.(sim.value[:, 1])
    μ̂ = logistic.(sim.value[:, 2])
    ν̂ = exp.(sim.value[:, 3])

    return(b̂, μ̂, ν̂)
end