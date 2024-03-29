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
"""
    location(pd::PearsonType1)

Obtain distribution location parameter a.
"""
function location(pd::PearsonType1)
    return pd.a    
end

"""
    scale(pd::PearsonType1)

Obtain distribution scale, given by b-a.
"""
function scale(pd::PearsonType1)
    return pd.b - pd.a
end

"""
    shape(pd::PearsonType1)

Obtain distribution shape parameters α and β.
"""
function shape(pd::PearsonType1)
    return (pd.α, pd.β)
end

"""
    shape(pd::PearsonType1)

Obtain all distribution parameters a, b, α and β.
"""
function params(pd::PearsonType1)
 return pd.a, pd.b, shape(pd)...
end



# Evaluations
"""
    cdf(pd::PearsonType1, x::Real)

Compute the cumulative distribution function value of pd at point x.
"""
function cdf(pd::PearsonType1, x::Real)
    return cdf(Beta(shape(pd)...),(x-location(pd))/scale(pd))
end

"""
    insupport(pd::PearsonType1, x::Real)

Establish if the point x is within the support of the distribution pd.
"""
function insupport(pd::PearsonType1, x::Real)
    return minimum(pd) <= x <= maximum(pd)
end

"""
    logpdf(pd::PearsonType1, x::Real)

Compute the log of the value of the probability density function of pd at point x.
"""
function logpdf(pd::PearsonType1, x::Real)
    return logpdf(Beta(shape(pd)...), (x-location(pd))/scale(pd)) - log(scale(pd))
end

"""
    maximum(pd::PearsonType1)

Obtain the upper limit of the distribution pd, b.
"""
function maximum(pd::PearsonType1)
    return params(pd)[2]
end

"""
    minimum(pd::PearsonType1)

Obtain the lower limit of the distribution pd, a.
"""
function minimum(pd::PearsonType1)
    return params(pd)[1]
end

"""
    quantile(pd::PearsonType1, p::Real)

Compute the quantile of probability p.
"""
function quantile(pd::PearsonType1, p::Real)
    return location(pd) + scale(pd) * quantile(Beta(shape(pd)...), p)
end

"""
    rand(rng::Random.AbstractRNG, pd::PearsonType1)

Generate a random realization of the distribution pd.
"""
function rand(rng::Random.AbstractRNG, pd::PearsonType1)
    return location(pd) + scale(pd) * rand(rng, Beta(shape(pd)...))
end



# Statistics
"""
    mean(pd::PearsonType1)

Obtain the expectation of the distribution pd.
"""
function mean(pd::PearsonType1)
    return location(pd) + scale(pd) * mean(Beta(shape(pd)...))
end

"""
    var(pd::PearsonType1)

Obtain the variance of the distribution pd.
"""
function var(pd::PearsonType1)
    return scale(pd)^2 * var(Beta(shape(pd)...))
end

"""
    std(pd::PearsonType1)

Obtain the standard deviation of the distribution pd.
"""
function std(pd::PearsonType1)
    return scale(pd) * std(Beta(shape(pd)...))
end

"""
    modes(pd::PearsonType1)

Obtain all modes (if this makes sense) of the distribution pd.
"""
function modes(pd::PearsonType1)
    return location(pd) .+ scale(pd) .* modes(Beta(shape(pd)...))
end

"""
    mode(pd::PearsonType1)

Obtain the first mode of the distribution pd.
"""
function mode(pd::PearsonType1)
    return location(pd) + scale(pd)* mode(Beta(shape(pd)...))
end

"""
    skewness(pd::PearsonType1)

Obtain the skewness of the distribution pd.
"""
function skewness(pd::PearsonType1)
    return skewness(Beta(shape(pd)...))
end

"""
    kurtosis(pd::PearsonType1)
    kurtosis(pd::PearsonType1, correction::Bool)

Obtain the excess kurtosis (if correction = false) or kurtosis (if correction = true) of the distribution pd .
"""
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

"""
    entropy(pd::PearsonType1)
    entropy(pd::PearsonType1, base::Real)

Compute the entropy value of distribution pd, w.r.t. a given base if given.
"""
function entropy(pd::PearsonType1)
    return entropy(Beta(shape(pd)...)) + log(scale(pd))
end

function entropy(pd::PearsonType1, base::Real)
    return entropy(pd)/log(base)
end



# MME
"""
    fit_mme(pd::Type{<:PearsonType1}, y::Vector{<:Real}, a::Real)
    fit_mme(pd::Type{<:PearsonType1}, y::Vector{<:Real})

Estimate parameters of a PearsonType1 distribution with method of moments.

The location parameter a can be fixed or not.
"""

function fit_mme(pd::Type{<:PearsonType1}, y::Vector{<:Real}, a::Real)
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
        @warn "The parameters associated with this data cannot be found by method of moments"
        return nothing
    end
    
    # parameters estimations
    sca = a2 - a1
    α = m1 + 1
    β = m2 + 1
    a = a
    b = a + sca

    # support verification
    if a >= minimum(y) 
        a = minimum(y) 
    end 
    if b <= maximum(y) 
        b = maximum(y) 
    end #-> message d'erreur

    return PearsonType1(a, b, α, β)
end

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
        @warn "The parameters associated with this data cannot be found by method of moments"
        return nothing
    end
    
    # parameters estimations
    sca = a2 - a1
    α = m1 + 1
    β = m2 + 1
    a = mm - sca * α/(α+β)
    b = a + sca

    # support verification
    if a >= minimum(y) 
        a = minimum(y) 
    end
    if b <= maximum(y) 
        b = maximum(y) 
    end #-> message d'erreur

    return PearsonType1(a, b, α, β)
end



# MLE
"""
    fit_mle(pd::Type{<:PearsonType1}, y::Vector{<:Real}, initialvalues::Vector{<:Real}, a::Real)
    fit_mle(pd::Type{<:PearsonType1}, y::Vector{<:Real}, initialvalues::Vector{<:Real})
    fit_mle(pd::Type{<:PearsonType1}, y::Vector{<:Real}, a::Real)
    fit_mle(pd::Type{<:PearsonType1}, y::Vector{<:Real})

Estimate parameters of a PearsonType1 distribution with likelihood maximization.

The location parameter a can be fixed or not. When a is given in argument, the initialvalues vector has to be of size 0 or 3. If a is not given, the initialvalues vector has to be of size 0 or 4
"""

function fit_mle(pd::Type{<:PearsonType1}, y::Vector{<:Real}, initialvalues::Vector{<:Real}, a::Real)
 
    iv = [initialvalues[1], initialvalues[2], initialvalues[3]]
    initialvalues[1] = max(initialvalues[1], maximum(y)) + .01*maximum(y)

    loglike(θ::Vector{<:Real}) = sum(logpdf.(Beta(θ[2], θ[3]), (y.-a)./(θ[1]-a)) .- log(θ[1]-a))
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
        
    return PearsonType1(a, θ̂[1], θ̂[2], θ̂[3])
end

function fit_mle(pd::Type{<:PearsonType1}, y::Vector{<:Real}, initialvalues::Vector{<:Real})
 
    iv = [initialvalues[1], initialvalues[2], initialvalues[3], initialvalues[4]]
    initialvalues[1] = min(initialvalues[1], minimum(y)) - abs(.01*minimum(y))
    initialvalues[2] = max(initialvalues[2], maximum(y)) + .01*maximum(y)

    loglike(θ::Vector{<:Real}) = sum(logpdf.(PearsonType1(θ...),y))
    fobj(θ) = -loglike(θ)

    lower = [-Inf, maximum(y), 2*eps(), 2*eps()]
    upper = [minimum(y), Inf, Inf, Inf]
    
    res = optimize(fobj, lower, upper, initialvalues, autodiff = :forward)
    
    if Optim.converged(res)
        θ̂ = Optim.minimizer(res)
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        θ̂ = iv
    end
        
    return PearsonType1(θ̂...)
end

function fit_mle(pd::Type{<:PearsonType1}, y::Vector{<:Real}, a::Real)
    initialvalues =  getinitialvalues(PearsonType1, y, a)
    fd = fit_mle(pd, y, initialvalues[2:4], a)
    return fd
end 

function fit_mle(pd::Type{<:PearsonType1}, y::Vector{<:Real})
    initialvalues =  getinitialvalues(PearsonType1, y)
    fd = fit_mle(pd, y, initialvalues)
    return fd
end 



# Mixed method
"""
    getinitialvalues(pd::Type{<:PearsonType1}, y::Vector{<:Real}, a::Real)
    getinitialvalues(pd::Type{<:PearsonType1}, y::Vector{<:Real})

Find initial values for the likelihood maximization algorithm. 

Initial values of shape parameters α and β are estimated by method of moments and initial values of boundary parameters a (when a is not given) and b are estimated by likelihood maximization.
"""

function getinitialvalues(pd::Type{<:PearsonType1}, y::Vector{<:Real}, a::Real)
    α, β = shape(fit_mme(pd, y))
    a = min(a, minimum(y)-abs(.01*minimum(y)))
    initialvalue = [maximum(y)+.01*maximum(y)]

    loglike(θ::Vector{<:Real}) = sum(logpdf.(Beta(α, β), (y.-a)./(θ[1]-a)) .- log(θ[1]-a))

    fobj(θ) = -loglike(θ)

    lower = [maximum(y)]
    upper = [Inf]

    res = optimize(fobj, lower, upper, initialvalue, autodiff = :forward)

    if Optim.converged(res)
        θ̂ = Optim.minimizer(res)
        return [a, θ̂[1], α, β]
    else
        @warn "The getinitialvalues algorithm did not find a solution. Maybe try with different initial values or with another method."
    end
end

function getinitialvalues(pd::Type{<:PearsonType1}, y::Vector{<:Real})
    α, β = shape(fit_mme(pd, y))
    initialvalues = [minimum(y)- abs(.01*minimum(y)), maximum(y)+.01*maximum(y)]

    loglike(θ::Vector{<:Real}) = sum(logpdf.(Beta(α, β), (y.-θ[1])./(θ[2]-θ[1])) .- log(θ[2]-θ[1]))

    fobj(θ) = -loglike(θ)

    lower = [-Inf, maximum(y)]
    upper = [minimum(y), Inf]

    res = optimize(fobj, lower, upper, initialvalues, autodiff = :forward)

    if Optim.converged(res)
        θ̂ = Optim.minimizer(res)
        return [θ̂[1], θ̂[2], α, β]
    else
        @warn "The getinitialvalues algorithm did not find a solution. Maybe try with different initial values or with another method."
    end
end



# Bayesian fitting
"""
    fit_bayes(pd::Type{<:PearsonType1}, y::Vector{<:Real}, prior::Real, niter::Int, warmup::Int, a::Real)

Estimate parameters of a PearsonType1b distribution with Bayesian inference.

Use NUTS (No U-Turn) sampler of MambaLite.
"""

function fit_bayes(pd::Type{<:PearsonType1}, y::Vector{<:Real}, prior::Real, niter::Int, warmup::Int, a::Real)
    nparam = 3
    iv = log.(getinitialvalues(PearsonType1, y))
    initialvalues = iv[2:4]

    # Defines the llh function and the gradient for the NUTS algo
    logf(θ::DenseVector) = sum(logpdf.(PearsonType1(a, exp(θ[1]), exp(θ[2]), exp(θ[3])), y)) + logpdf(Exponential(prior), exp(θ[1]))
    Δlogf(θ::DenseVector) = ForwardDiff.gradient(logf, θ)
    function logfgrad(θ::DenseVector)
        ll = logf(θ)
        g = Δlogf(θ)
        return ll, g
    end

    # MCMC
    sim = Chains(niter, nparam, start=(warmup+1))
    θ = NUTSVariate(initialvalues, logfgrad)
    for i in 1:niter
        MambaLite.sample!(θ, adapt=(i<=warmup))
        if i>warmup
            sim[i, :, 1] = θ
        end
    end

    b̂ = exp.(sim.value[:, 1])
    α̂ = exp.(sim.value[:, 2])
    β̂ = exp.(sim.value[:, 3])

    return(b̂, α̂, β̂)
end

#function fit_bayes(pd::Type{<:PearsonType1}, y::Vector{<:Real}, prior::Real, niter::Int, warmup::Int)
#    nparam = 4
#    iv = getinitialvalues(PearsonType1, y)
#    initialvalues = [iv[1], log(iv[2]), log(iv[3]), log(iv[4])]
#
#    # Defines the llh function and the gradient for the NUTS algo
#    logf(θ::DenseVector) = sum(logpdf(PearsonType1(θ[1], exp(θ[2]), exp(θ[3]), exp(θ[4])), y)) + logpdf(Exponential(prior), exp(θ[2]))
#    Δlogf(θ::DenseVector) = ForwardDiff.gradient(logf, θ)
#    function logfgrad(θ::DenseVector)
#        ll = logf(θ)
#        g = Δlogf(θ)
#        return ll, g
#    end
#
#    # MCMC
#    sim = Chains(niter, nparam, start=(warmup+1))
#    θ = NUTSVariate(initialvalues, logfgrad)
#    for i in 1:niter
#        MambaLite.sample!(θ, adapt=(i<=warmup))
#        if i>warmup
#            sim[i, :, 1] = θ
#        end
#    end
#
#    â = sim.value[:, 1]
#    b̂ = exp.(sim.value[:, 2])
#    α̂ = exp.(sim.value[:, 3])
#    β̂ = exp.(sim.value[:, 4])
#    return(â, b̂, α̂, β̂)
#end