"""
    PearsonType1(a,b,α,β)

The *Pearson Type 1 distribution* with shape parameters `α` and `β` defined on the interval (`a`, `b`) has the probability density function for ``a<y<b``
```math
f(y; a, b, \\alpha, \\beta) = \\frac{1}{B(\\alpha, \\beta)} \\frac{(y-a)^{\\alpha-1} (b-y)^{\\beta-1}}{(b-a)^{\\alpha+\\beta-1}},
```
and 0 elsewhere.


```julia
PearsonType1()   # Pearson Type 1 distribution on the unit interval with shape parameters (1,1) i.e. Uniform(0, 1)


params(d)        # Get the parameters, i.e. (a, b, α, β)

location(d)      # Get the location parameter, i.e. a
scale(d)         # Get the scale parameter, i.e. b-a
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

PearsonType1(a::Real, b::Real, α::Real, β::Real; check_args::Bool=true) = PearsonType1(promote(a, b, α , β)...; check_args=check_args)
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
    return cdf(Beta(shape(pd)...),(x-location(pd))/scale(pd))
end

function insupport(pd::PearsonType1, x::Real)
    return minimum(pd) <= x <= maximum(pd)
end

function logpdf(pd::PearsonType1, x::Real)
    return logpdf(Beta(shape(pd)...), (x-location(pd))/scale(pd)) - log(scale(pd))
end

function maximum(pd::PearsonType1)
    return params(pd)[2]
end

function minimum(pd::PearsonType1)
    return params(pd)[1]
end

function quantile(pd::PearsonType1, p::Real)
    return location(pd) + scale(pd) * quantile(Beta(shape(pd)...), p)
end

function rand(rng::Random.AbstractRNG, pd::PearsonType1)
    return location(pd) + scale(pd) * rand(rng, Beta(shape(pd)...))
end



# Statistics

function mean(pd::PearsonType1)
    return location(pd) + scale(pd) * mean(Beta(shape(pd)...))
end

function var(pd::PearsonType1)
    return scale(pd)^2 * var(Beta(shape(pd)...))
end

function std(pd::PearsonType1)
    return scale(pd) * std(Beta(shape(pd)...))
end

function modes(pd::PearsonType1)
    return location(pd) .+ scale(pd) .* modes(Beta(shape(pd)...))
end

function mode(pd::PearsonType1)
    return location(pd) + scale(pd)* mode(Beta(shape(pd)...))
end

function skewness(pd::PearsonType1)
    return skewness(Beta(shape(pd)...))
end

function kurtosis(pd::PearsonType1) # excess kurtosis
    return kurtosis(Beta(shape(pd)...))
end

function kurtosis(pd::PearsonType1, correction::Bool) # kurtosis
    td = Beta(shape(pd)...)
    if correction
        return kurtosis(td)
    end
    return kurtosis(td) + 3.0
end

function entropy(pd::PearsonType1)
    return entropy(Beta(shape(pd)...)) + log(scale(pd))
end

function entropy(pd::PearsonType1, base::Real)
    return entropy(pd)/log(base)
end



# fit by method of moments

function fit_mme(pd::Type{<:PearsonType1}, y::Vector{<:Real})
    # sample moments
    mm = mean(y)
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
    sca = a2 - a1
    α = m1 + 1
    β = m2 + 1
    a = mm - sca * α/(α+β)
    b = a + sca

    # vérifier le support
    if a > minimum(y) 
        a = minimum(y) 
    end #-> message d'erreur 
    if b < maximum(y) 
        b = maximum(y) 
    end #-> message d'erreur

    return PearsonType1(a, b, α, β)
end



# fit by maximum likelihood 

function fit_mle(pd::Type{<:PearsonType1}, y::Vector{<:Real}, initialvalues::Vector{<:Real})
 
#     # PearsonDS propose de +-0.1 aux valeurs initiales
#     if initialvalues[2] > 0
#         initialvalues[1] = min(initialvalues[1], minimum(y)) # - 0.1
#         initialvalues[2] = max(initialvalues[2], maximum(y)) # + 0.1
#     else
#         initialvalues[1] = max(initialvalues[1], maximum(y)) # + 0.1
#         initialvalues[2] = min(initialvalues[2], minimum(y)) # - 0.1
#     end 

    loglike(θ::Vector{<:Real}) = -sum(logpdf.(PearsonType1(θ...),y))
    fobj(θ) = -loglike(θ)

    lower = [-Inf, maximum(y), 0, 0]
    upper = [minimum(y), Inf, Inf, Inf]
    
    res = optimize(fobj, lower, upper, initialvalues, autodiff = :forward)
    
    if Optim.converged(res)
        θ̂ = Optim.minimizer(res)
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        θ̂ = initialvalues
    end
        
    return PearsonType1(θ̂...)
end

function fit_mle(pd::Type{<:PearsonType1}, y::Vector{<:Real})
    fd = PMP.fit_mle(pd, y, [params(fit_mme(PearsonType1, y))...])
    return fd
end 



# find initial values for fit_mle

function getinitialvalues(pd::Type{<:PearsonType1}, y::Vector{<:Real})
    α, β = shape(fit_mme(pd, y))
    a = minimum(y)
    b = maximum(y)

    loglike(θ₁::Real, θ₂::Real) = sum(logpdf.(PearsonType1(θ₁, θ₂, α, β),y))
    fobj(θ₁, θ₂) = -loglike(θ₁, θ₂)

    lower = [-Inf, maximum(y)]
    upper = [minimum(y), Inf]

    res = optimize(fobj, lower, upper, [a, b], autodiff = :forward)

    if Optim.converged(res)
        θ̂₁, θ̂₂ = Optim.minimizer(res)
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        θ̂₁, θ̂₂ = initialvalues
    end

    return θ̂₁, θ̂₂, α, β
end

# trouver alpha beta par moment, a b par mle -> getinitialvalues pour ensuite faire mle pour chaque params