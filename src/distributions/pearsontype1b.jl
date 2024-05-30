"""
    PearsonType1b(b,α,β)

The *Pearson Type 1 distribution* with shape parameters `α` and `β` defined on the interval (0, `b`) has the probability density function for ``0<y<b``
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
    @check_args PearsonType1b (b, zero(b) < b) (α, zero(α) < α) (β, zero(β) < β)
    return PearsonType1b{T}(b, α, β)
end

PearsonType1b(b::Real, α::Real, β::Real; check_args::Bool=true) = PearsonType1b(promote(b, α, β)...; check_args=check_args)
PearsonType1b(b::Real, α::Integer, β::Integer; check_args::Bool=true) = PearsonType1b(float(b), float(α), float(β); check_args=check_args)
PearsonType1b() = PearsonType1b{Float64}(1.0, 1.0, 1.0)



# Parameters
"""
    location(pd::PearsonType1b)

Obtain distribution location parameter 0.
"""
function location(pd::PearsonType1b)
    return 0    
end

"""
    scale(pd::PearsonType1b)

Obtain distribution scale parameter b.
"""
function scale(pd::PearsonType1b)
    return pd.b
end

"""
    shape(pd::PearsonType1b)

Obtain distribution shape parameters α and β.
"""
function shape(pd::PearsonType1b)
    return (pd.α, pd.β)
end

"""
    shape(pd::PearsonType1b)

Obtain all distribution parameters b, α and β.
"""
function params(pd::PearsonType1b)
 return pd.b, shape(pd)...
end



# Evaluations
"""
    cdf(pd::PearsonType1b, x::Real)

Compute the cumulative distribution function value of pd at point x.
"""
function cdf(pd::PearsonType1b, x::Real)
    return cdf(Beta(shape(pd)...), x/pd.b)
end

"""
    insupport(pd::PearsonType1b, x::Real)

Establish if the point x is within the support of the distribution pd.
"""
function insupport(pd::PearsonType1b, x::Real)
    return 0 <= x <= maximum(pd)
end

"""
    logpdf(pd::PearsonType1b, x::Real)

Compute the log of the value of the probability density function of pd at point x.
"""
function logpdf(pd::PearsonType1b, x::Real)
    return logpdf(Beta(shape(pd)...), x/pd.b) - log(pd.b)
end

"""
    maximum(pd::PearsonType1b)

Obtain the upper limit of the distribution pd, b.
"""
function maximum(pd::PearsonType1b)
    return scale(pd)
end

"""
    minimum(pd::PearsonType1b)

Obtain the lower limit of the distribution pd, 0.
"""
function minimum(pd::PearsonType1b)
    return 0
end

"""
    quantile(pd::PearsonType1b, p::Real)

Compute the quantile of probability p.
"""
function quantile(pd::PearsonType1b, p::Real)
    return pd.b * quantile(Beta(shape(pd)...), p)
end

"""
    rand(rng::Random.AbstractRNG, pd::PearsonType1b)

Generate a random realization of distribution pd.
"""
function rand(rng::Random.AbstractRNG, pd::PearsonType1b)
    return pd.b * rand(rng, Beta(shape(pd)...))
end



# Statistics
"""
    mean(pd::PearsonType1b)

Obtain the expectation of the distribution pd.
"""
function mean(pd::PearsonType1b)
    return pd.b * mean(Beta(shape(pd)...))
end

"""
    var(pd::PearsonType1b)

Obtain the variance of the distribution pd.
"""
function var(pd::PearsonType1b)
    return pd.b^2 * var(Beta(shape(pd)...))
end

"""
    std(pd::PearsonType1b)

Obtain the standard deviation of the distribution pd.
"""
function std(pd::PearsonType1b)
    return pd.b * std(Beta(shape(pd)...))
end

"""
    modes(pd::PearsonType1b)

Obtain all modes (if this makes sense) of the distribution pd.
"""
function modes(pd::PearsonType1b)
    return pd.b .* modes(Beta(shape(pd)...))
end

"""
    mode(pd::PearsonType1b)

Obtain the first mode of the distribution pd.
"""
function mode(pd::PearsonType1b)
    return pd.b * mode(Beta(shape(pd)...))
end

"""
    skewness(pd::PearsonType1b)

Obtain the skewness of the distribution pd.
"""
function skewness(pd::PearsonType1b)
    return skewness(Beta(shape(pd)...))
end

"""
    kurtosis(pd::PearsonType1b)
    kurtosis(pd::PearsonType1b, correction::Bool)

Obtain the excess kurtosis (if correction = false) or kurtosis (if correction = true) of the distribution pd .
"""
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

"""
    entropy(pd::PearsonType1b)
    entropy(pd::PearsonType1b, base::Real)

Compute the entropy value of distribution pd, w.r.t. a given base if given.
"""
function entropy(pd::PearsonType1b)
    return entropy(Beta(shape(pd)...)) + log(pd.b)
end

function entropy(pd::PearsonType1b, base::Real)
    return entropy(pd)/log(base)
end



# MME
"""
    fit_mme(pd::Type{<:PearsonType1b}, y::Vector{<:Real})

Estimate parameters of a PearsonType1b distribution with method of moments.
"""

function fit_mme(pd::Type{<:PearsonType1b}, y::Vector{<:Real})
    # sample moment
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



# MLE
"""
    fit_mle(pd::Type{<:PearsonType1b}, y::Vector{<:Real}, initialvalues::Vector{<:Real})
    fit_mle(pd::Type{<:PearsonType1b}, y::Vector{<:Real})

Estimate parameters of a PearsonType1b distribution with likelihood maximization.
"""

function fit_mle(pd::Type{<:PearsonType1b}, y::Vector{<:Real}, initialvalues::Vector{<:Real})
 
    iv = [initialvalues[1], initialvalues[2], initialvalues[3]]
    initialvalues[1] = max(initialvalues[1], maximum(y)) + .01*maximum(y)

    loglike(θ::Vector{<:Real}) = sum(logpdf.(PearsonType1b(θ...),y))
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
        
    return PearsonType1b(θ̂...)
end

function fit_mle(pd::Type{<:PearsonType1b}, y::Vector{<:Real})
    initialvalues =  getinitialvalues(PearsonType1b, y)
    fd = PMP.fit_mle(pd, y, initialvalues)
    return fd
end 



# Mixed method
"""
    getinitialvalues(pd::Type{<:PearsonType1b}, y::Vector{<:Real})

Find initial values for the likelihood maximization algorithm. 

Initial values of shape parameters α and β are estimated by method of moments and initial value of scale parameter b is estimated by likelihood maximization.
"""

function getinitialvalues(pd::Type{<:PearsonType1b}, y::Vector{<:Real})
    α, β = shape(fit_mme(pd, y))
    initialvalue = [maximum(y) + .01*maximum(y)]

    loglike(θ::Vector{<:Real}) = sum(logpdf.(Beta(α, β), y./θ[1]) .- log(θ[1]))
    fobj(θ) = -loglike(θ)

    lower = [maximum(y)]
    upper = [Inf]

    res = optimize(fobj, lower, upper, initialvalue, autodiff = :forward)

    if Optim.converged(res)
        θ̂ = Optim.minimizer(res)
    else
        @warn "The getinitialvalues algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        θ̂ = initialvalue
    end

    return [θ̂[1], α, β]
end



# Bayesian fitting
"""
    fit_bayes(pd::Type{<:PearsonType1b}, y::Vector{<:Real}, prior::Real, niter::Int, warmup::Int)

Estimate parameters of a PearsonType1b distribution with Bayesian inference.

Use NUTS (No U-Turn) sampler. The prior refers to the prior distribution of the bound parameter b.
"""

function fit_bayes(pd::Type{<:PearsonType1b}, prior::ContinuousUnivariateDistribution, y::Vector{<:Real}, niter::Int, warmup::Int)
    nparam = 3
    initialvalues = log.(getinitialvalues(pd, y))

    # Defines the llh function and the gradient for the NUTS algo
    logf(θ::DenseVector) = sum(logpdf.(pd(exp(θ[1]), exp(θ[2]), exp(θ[3])), y)) + logpdf(prior, exp(θ[1]))
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
    α̂ = exp.(sim.value[:, 2])
    β̂ = exp.(sim.value[:, 3])

    return(b̂, α̂, β̂)
end



"""
    fit_bayes_MH(pd::Type{<:PearsonType1b}, y::Vector{<:Real}, π_b::ContinuousUnivariateDistribution, π_α::ContinuousUnivariateDistribution, π_β::ContinuousUnivariateDistribution, warmup::Int=5000, thin::Int=10, niter::Int=20000)

Estimate parameters of a PearsonType1b distribution with Metropolis-Hasting algorithm.
"""

function fit_bayes_MH(pd::Type{<:PearsonType1b}, y::Vector{<:Real}, 
    π_b::ContinuousUnivariateDistribution, 
    π_α::ContinuousUnivariateDistribution, 
    π_β::ContinuousUnivariateDistribution,
    warmup::Int=5000, 
    thin::Int=10, 
    niter::Int=20000)

    b̂ = Float64[]
    α̂ = Float64[]
    β̂ = Float64[]

    initialvalues =  getinitialvalues(PearsonType1b, y)
    σ = initialvalues./20

    push!(b̂, initialvalues[1])
    push!(α̂, initialvalues[2])
    push!(β̂, initialvalues[3])

    acc = falses(3, niter)

    for it in 1:niter
            # for b
    
        q_b = rand(Normal(b̂[end],σ[1]))
        log_r = - sum(logpdf(PearsonType1b(b̂[end], α̂[end], β̂[end]), y)) - logpdf(π_b, b̂[end]) + sum(logpdf(PearsonType1b(q_b, α̂[end], β̂[end]), y)) + logpdf(π_b, q_b)
        cond = log(rand(Uniform(0,1))) < log_r 
        acc[1, it] = cond
        b_temp = q_b*cond + b̂[end]*(1-cond)
        b̂ = push!(b̂, b_temp)
        
        # for α
        q_α = rand(Normal(α̂[end],σ[2]))
        log_r = - sum(logpdf(PearsonType1b(b̂[end], α̂[end], β̂[end]), y)) - logpdf(π_α, α̂[end]) + sum(logpdf(PearsonType1b(b̂[end], q_α, β̂[end]), y)) + logpdf(π_α, q_α)
        cond = log(rand(Uniform(0,1))) < log_r
        acc[2, it] = cond
        α_temp = q_α*cond + α̂[end]*(1-cond)
        α̂ = push!(α̂, α_temp)
    
        # for β
        q_β = rand(Normal(β̂[end],σ[3]))
        log_r = - sum(logpdf(PearsonType1b(b̂[end], α̂[end], β̂[end]), y)) - logpdf(π_β, β̂[end]) + sum(logpdf(PearsonType1b(b̂[end], α̂[end], q_β), y)) + logpdf(π_β, q_β)
        cond = log(rand(Uniform(0,1))) < log_r
        acc[3, it] = cond
        β_temp = q_β*cond + β̂[end]*(1-cond)
        β̂ = push!(β̂, β_temp)
    
        # updating instrumental distribution
        if it % 50 == 0
            if (it<=warmup)
                accrate = vec(mean(acc[:, it-50+1:it], dims=2))
                σ = update_stepsize.(σ, accrate)
            end
        end
    end

    b̂ = b̂[warmup:thin:niter]
    α̂ = α̂[warmup:thin:niter]
    β̂ = β̂[warmup:thin:niter]

    return(b̂, α̂, β̂)
end