"""
    PearsonType1b(b,α,β)

The *Pearson Type 1 distribution* with shape parameters `α` and `β` defined on the interval (0, `b`) has the probability density function for ``a<y<b``
```math
f(y; b, \\alpha, \\beta) = \\frac{1}{B(\\alpha, \\beta)} \\frac{y^{\\alpha-1} (b-y)^{\\beta-1}}{b^{\\alpha+\\beta-1}},
```
and 0 elsewhere.


```julia
PearsonType1b()   # Pearson Type 1 distribution on the unit interval with shape parameters (1,1) i.e. Uniform(0, 1)


params(d)        # Get the parameters, i.e. (b, α, β)

scale(d)         # Get the scale parameter, i.e. b
shape(d)         # Get the shape parameters, i.e. (α, β)
```

External links

* [Pearson Type 1 distribution on Wikipedia](https://en.wikipedia.org/wiki/Pearson_distribution#The_Pearson_type_I_distribution)

"""
struct PearsonType1b{T<:Real} <: ContinuousUnivariateDistribution
    b::T
    α::T
    β::T
    PearsonType1b{T}(b::T, α::T, β::T) where {T<:Real} = new{T}(b, α, β)
end

function PearsonType1b(b::T, α::T, β::T ; check_args::Bool=true) where {T <: Real}
    @check_args PearsonType1b (b, b>0) (α, α > zero(α)) (β, β > zero(β))
    return PearsonType1b{T}(b, α, β)
end

PearsonType1b(b::Real, α::Real, β::Real; check_args::Bool=true) = PearsonType1b(promote(b, α , β)...; check_args=check_args)
PearsonType1b(b::Real, α::Integer, β::Integer; check_args::Bool=true) = PearsonType1b(float(b), float(α), float(β); check_args=check_args)
PearsonType1b() = PearsonType1b{Float64}(1.0, 1.0, 1.0)



# Parameters

function location(pd::PearsonType1b)
    return 0    
end

function scale(pd::PearsonType1b)
    return pd.b
end

function shape(pd::PearsonType1b)
    return (pd.α, pd.β)
end

function params(pd::PearsonType1b)
 return pd.b, shape(pd)...
end



# Evaluations

function cdf(pd::PearsonType1b, x::Real)
    return cdf(Beta(shape(pd)...), x/pd.b)
end

function insupport(pd::PearsonType1b, x::Real)
    return 0 <= x <= maximum(pd)
end

function logpdf(pd::PearsonType1b, x::Real)
    return logpdf(Beta(shape(pd)...), x/pd.b) - log(pd.b)
end

function maximum(pd::PearsonType1b)
    return scale(pd)
end

function minimum(pd::PearsonType1b)
    return 0
end

function quantile(pd::PearsonType1b, p::Real)
    return pd.b * quantile(Beta(shape(pd)...), p)
end

function rand(rng::Random.AbstractRNG, pd::PearsonType1b)
    return pd.b * rand(rng, Beta(shape(pd)...))
end



# Statistics

function mean(pd::PearsonType1b)
    return pd.b * mean(Beta(shape(pd)...))
end

function var(pd::PearsonType1b)
    return pd.b^2 * var(Beta(shape(pd)...))
end

function std(pd::PearsonType1b)
    return pd.b * std(Beta(shape(pd)...))
end

function modes(pd::PearsonType1b)
    return pd.b .* modes(Beta(shape(pd)...))
end

function mode(pd::PearsonType1b)
    return pd.b * mode(Beta(shape(pd)...))
end

function skewness(pd::PearsonType1b)
    return skewness(Beta(shape(pd)...))
end

function kurtosis(pd::PearsonType1b) # excess kurtosis
    return kurtosis(Beta(shape(pd)...))
end

function kurtosis(pd::PearsonType1b, correction::Bool) # kurtosis
    td = Beta(shape(pd)...)
    if correction
        return kurtosis(td)
    end
    return kurtosis(td) + 3.0
end

function entropy(pd::PearsonType1b)
    return entropy(Beta(shape(pd)...)) + log(pd.b)
end

function entropy(pd::PearsonType1b, base::Real)
    return entropy(pd)/log(base)
end



# fit by method of moments
function fit_mme(pd::Type{<:PearsonType1b}, y::Vector{<:Real})
    # sample moments
    vv = var(y)
    ss = skewness(y)
    kk = kurtosis(y) + 3 # not excess kurtosis

    # the kurtosis is bounded below by the squared skewness plus 1
    if ss^2 > kk-1 
        @error "There are no probability distributions with these moments" 
    end

    aa = 2*kk - 3*ss^2 - 6
    bb = ss*(kk + 3)
    cc = 4*kk - 3*ss^2

    a1 = sqrt(vv)/2 * ((-bb-sqrt(bb^2-4*cc*aa))/aa)
    a2 = sqrt(vv)/2 * ((-bb+sqrt(bb^2-4*cc*aa))/aa)
    if a1 > 0
        tmp = a1 
        a1 = a2 
        a2 = tmp
    end

    m1 = -(bb+a1*(10*kk-12*ss^2-18)/sqrt(vv)) / (sqrt(bb^2-4*cc*aa))
    m2 = -(-bb-a2*(10*kk-12*ss^2-18)/sqrt(vv)) / (sqrt(bb^2-4*cc*aa))
    #@assert m1 > -1 && m2 > -1
    if m1 < -1 || m2 < -1
        return nothing
    end
    
    # parameters estimations
    b = a2 - a1
    α = m1 + 1
    β = m2 + 1

    return PearsonType1b(b, α, β)
end



# fit by maximum likelihood 

function fit_mle(pd::Type{<:PearsonType1b}, y::Vector{<:Real}, initialvalues::Vector{<:Real})
 
    # PearsonDS propose de +-0.1 aux valeurs initiales
    initialvalues[1] = max(initialvalues[2], maximum(y)) # + 0.1

    loglike(θ::Vector{<:Real}) = sum(logpdf.(PearsonType1b(θ...),y))

    fobj(θ) = -loglike(θ)

    res = optimize(fobj, initialvalues)

    if Optim.converged(res)
        θ̂ = Optim.minimizer(res)
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        θ̂ = Optim.minimizer(res)
    end
    
    return PearsonType1b(θ̂...)
end

# ne fonctionne pas : erreur 
# MethodError: no method matching setindex!(::NTuple{4, Float64}, ::Float64, ::Int64)
function fit_mle(pd::Type{<:PearsonType1b}, y::Vector{<:Real})
    
    initialvalues = PMP.fit_mme(PearsonType1b, y)
    
    return PMP.fit_mle(pd, y, initialvalues)  # retour un PearsonType1 b
end
