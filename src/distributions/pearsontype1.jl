"""
    PearsonType1(a,b,α,β)

The *Pearson Type 1 distribution* with shape parameters `α` and `β` defined on the interval (`a`, `b`)has the probability density
function for ``a<y<b``

```math
f(y; a, b, \\alpha, \\beta) = \\frac{1}{\\Beta(\\alpha, \\beta)} \\frac{(y-a)^{\\alpha-1} (b-y)^{\\beta-1}}{(b-a)^{\\alpha+\\beta-1},
```
and 0 elsewhere.


```julia
PearsonType1()   # Pearson Type 1 distribution on the unit interval with shape parameters (1,1) i.e. Uniform(0, 1)


params(d)        # Get the parameters, i.e. (a, b, α, β)

location(d)      # Get the location parameter, i.e. a
scale(d)         # Get the scape parameter, i.e. b-a
shape(d)         # Get the shape parameters, i.e. (α, β)
```

External links

* [Pearson Type 1 distribution on Wikipedia](https://en.wikipedia.org/wiki/Pearson_distribution#The_Pearson_type_I_distribution)

"""
struct PearsonType1{T<:Real} <: ContinuousUnivariateDistribution
    a::T
    b::T
    α::T
    β::T
    PearsonType1{T}(a::T, b::T, α::T, β::T) where {T<:Real} = new{T}(a, b, α, β)
end

function PearsonType1(a::T, b::T, α::T, β::T ; check_args::Bool=true) where {T <: Real}
    @check_args PearsonType1 ((a,b), b>a) (α, α > zero(α)) (β, β > zero(β))
    return PearsonType1{T}(a, b, α, β )
end

PearsonType1(a::Real, b::Real, α::Real, β::Real; check_args::Bool=true) = PearsonType1(promote(a, b, α, β)...; check_args=check_args)
PearsonType1(a::Real, b::Real, α::Integer, β::Integer; check_args::Bool=true) = PearsonType1(float(a), float(b), float(α), float(β); check_args=check_args)
PearsonType1() = PearsonType1{Float64}(0.0, 1.0, 1.0, 1.0)



# Parameters

function location(pd::PearsonType1)
    return pd.a    
end

function scale(pd::PearsonType1)
    return pd.b - pd.a
end

function shape(pd::PearsonType1)
    return (pd.α, pd.β)
end

function params(pd::PearsonType1)
 return pd.a, pd.b, shape(pd)...
end



# Evaluations

function cdf(pd::PearsonType1, x::Real)   
    td = getdistribution(pd)
    return cdf(td, x)
end

function getdistribution(pd::PearsonType1)
   return LocationScale(location(pd), scale(pd), Beta(shape(pd)...)) 
end

function insupport(pd::PearsonType1, x::Real)
    return minimum(pd) <= x <= maximum(pd)
end

function logpdf(pd::PearsonType1, x::Real)
    td = getdistribution(pd)
    return logpdf(td, x)
end

function maximum(pd::PearsonType1)
    return params(pd)[2]
end

function minimum(pd::PearsonType1)
    return params(pd)[1]
end

function quantile(pd::PearsonType1, p::Real)
    td = getdistribution(pd)
    return quantile(td, p)
end

function rand(rng::Random.AbstractRNG, pd::PearsonType1)
    td = getdistribution(pd)
    return rand(rng, td)
end

# fit by maximum likelihood 

function fit_mle(pd::Type{<:PearsonType1}, y::Vector{<:Real}, initialvalues::Vector{<:Real})
    
    loglike(θ::Vector{<:Real}) = sum(logpdf.(PearsonType1(θ...),y))

    fobj(θ) = -loglike(θ)

    initialvalues = [-1., 1., 2., 3.]

    res = optimize(fobj, initialvalues)

    if Optim.converged(res)
        θ̂ = Optim.minimizer(res)
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        θ̂ = Optim.minimizer(res)
    end
    
    return PearsonType1(θ̂...)
    
end

function fit_mle(pd::Type{<:PearsonType1}, y::Vector{<:Real})

    # TODO Replace initial values with the estimations obtained with th method of moments.
    initialvalues = [minimum(y), maximum(y), 1., 1.]
    
    return fit_mle(pd, y, initialvalues)
    
end