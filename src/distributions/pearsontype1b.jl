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

    if ss>0 && m1>m2
        tmp = m1
        m1 = m2
        m2 = tmp
    end
    
    if ss<0 && m1<m2
        tmp = m1
        m1 = m2
        m2 = tmp
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
    fit_bayes(pd::Type{<:PearsonType1b}, y::Vector{<:Real}, π_b::ContinuousUnivariateDistribution; warmup::Int, niter::Int)
    fit_bayes(pd::Type{<:PearsonType1b}, y::Vector{<:Real}; warmup::Int, niter::Int)

Estimate parameters of a PearsonType1b distribution with Bayesian inference.

Use NUTS (No U-Turn) sampler. The prior refers to the prior distribution of the bound parameter b.
"""

function fit_bayes(pd::Type{<:PearsonType1b}, y::Vector{<:Real}, π_b::ContinuousUnivariateDistribution; warmup::Int, niter::Int)
    nparam = 3
    initialvalues = log.(getinitialvalues(pd, y))

    # Defines the llh function and the gradient for the NUTS algo
    logf(θ::DenseVector) = sum(logpdf.(pd(exp(θ[1]), exp(θ[2]), exp(θ[3])), y)) + logpdf(π_b, exp(θ[1]))
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
    
    b = exp.(sim.value[:, 1])
    α = exp.(sim.value[:, 2])
    β = exp.(sim.value[:, 3])

    return(b, α, β)
end

function fit_bayes(pd::Type{<:PearsonType1b}, y::Vector{<:Real}; warmup::Int, niter::Int)
    nparam = 3
    initialvalues = log.(getinitialvalues(pd, y))

    # Defines the llh function and the gradient for the NUTS algo
    logf(θ::DenseVector) = sum(logpdf.(pd(exp(θ[1]), exp(θ[2]), exp(θ[3])), y))
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
    
    b = exp.(sim.value[:, 1])
    α = exp.(sim.value[:, 2])
    β = exp.(sim.value[:, 3])

    return(b, α, β)
end



"""
    fit_bayes_MH(pd::Type{<:PearsonType1b}, y::Vector{<:Real}, π_b::ContinuousUnivariateDistribution, π_α::ContinuousUnivariateDistribution, π_β::ContinuousUnivariateDistribution; warmup::Int=5000, thin::Int=10, niter::Int=20000)
    fit_bayes_MH(pd::Type{<:PearsonType1b}, y::Vector{<:Real}; warmup::Int=5000, thin::Int=10, niter::Int=20000)

Estimate parameters of a PearsonType1b distribution with Metropolis-Hasting algorithm.
"""

function fit_bayes_MH(pd::Type{<:PearsonType1b}, y::Vector{<:Real}, 
    π_b::ContinuousUnivariateDistribution, 
    π_α::ContinuousUnivariateDistribution, 
    π_β::ContinuousUnivariateDistribution;
    warmup::Int=5000, thin::Int=10, niter::Int=20000)

    initialvalues = params(fit_mme(pd, y))
    δ = initialvalues./10

    acc = falses(3, niter)
    b = Vector{Float64}(undef, niter)
    α = Vector{Float64}(undef, niter)
    β = Vector{Float64}(undef, niter)
    
    b[1], α[1], β[1] = initialvalues
    
    for iter in 2:niter
        b̃ = b[iter-1] + randn()*δ[1]
        log_r = - loglikelihood(PearsonType1b(b[iter-1], α[iter-1], β[iter-1]), y) +
            loglikelihood(PearsonType1b(b̃, α[iter-1], β[iter-1]), y) +
            logpdf(π_b, b̃) - logpdf(π_b, b[iter-1])
        if log_r > log(rand())
            b[iter] = b̃
            acc[1, iter] = true
        else
            b[iter] = b[iter-1]
        end
        
        α̃ = exp(log(α[iter-1]) + randn()*δ[2])
        log_r = - loglikelihood(PearsonType1b(b[iter], α[iter-1], β[iter-1]), y) +
            loglikelihood(PearsonType1b(b[iter], α̃, β[iter-1]), y) +
            logpdf(π_α, α̃) - logpdf(π_α, α[iter-1])
        if log_r > log(rand())
            α[iter] = α̃
            acc[2, iter] = true
        else
            α[iter] = α[iter-1]
        end
        
        β̃ = exp(log(β[iter-1]) + randn()*δ[3])
        log_r = - loglikelihood(PearsonType1b(b[iter], α[iter], β[iter-1]), y) +
            loglikelihood(PearsonType1b(b[iter], α[iter], β̃), y) +
            logpdf(π_β, β̃) - logpdf(π_β, β[iter-1])
        if log_r > log(rand())
            β[iter] = β̃
            acc[3, iter] = true
        else
            β[iter] = β[iter-1]
        end
        
        # updating instrumental distribution
        if iter % 50 == 0
            if (iter<=warmup)
                accrate = vec(mean(acc[:, iter-50+1:iter], dims=2))
                δ = update_stepsize.(δ, accrate)
            end
        end
        
    end
    
    b = b[warmup:thin:niter]
    α = α[warmup:thin:niter]
    β = β[warmup:thin:niter]
    
    return b, α, β
end

function fit_bayes_MH(pd::Type{<:PearsonType1b}, y::Vector{<:Real}; warmup::Int=5000, thin::Int=10, niter::Int=20000)

    initialvalues = params(fit_mme(pd, y))
    δ = initialvalues./10

    acc = falses(3, niter)
    b = Vector{Float64}(undef, niter)
    α = Vector{Float64}(undef, niter)
    β = Vector{Float64}(undef, niter)
    
    b[1], α[1], β[1] = initialvalues
    
    for iter in 2:niter
        b̃ = b[iter-1] + randn()*δ[1]
        log_r = - loglikelihood(PearsonType1b(b[iter-1], α[iter-1], β[iter-1]), y) +
            loglikelihood(PearsonType1b(b̃, α[iter-1], β[iter-1]), y)
        if log_r > log(rand())
            b[iter] = b̃
            acc[1, iter] = true
        else
            b[iter] = b[iter-1]
        end
        
        α̃ = exp(log(α[iter-1]) + randn()*δ[2])
        log_r = - loglikelihood(PearsonType1b(b[iter], α[iter-1], β[iter-1]), y) +
            loglikelihood(PearsonType1b(b[iter], α̃, β[iter-1]), y)
        if log_r > log(rand())
            α[iter] = α̃
            acc[2, iter] = true
        else
            α[iter] = α[iter-1]
        end
        
        β̃ = exp(log(β[iter-1]) + randn()*δ[3])
        log_r = - loglikelihood(PearsonType1b(b[iter], α[iter], β[iter-1]), y) +
            loglikelihood(PearsonType1b(b[iter], α[iter], β̃), y)
        if log_r > log(rand())
            β[iter] = β̃
            acc[3, iter] = true
        else
            β[iter] = β[iter-1]
        end
        
        # updating instrumental distribution
        if iter % 50 == 0
            if (iter<=warmup)
                accrate = vec(mean(acc[:, iter-50+1:iter], dims=2))
                δ = update_stepsize.(δ, accrate)
            end
        end
        
    end
    
    b = b[warmup:thin:niter]
    α = α[warmup:thin:niter]
    β = β[warmup:thin:niter]
    
    return b, α, β
end



"""
    update_stepsize(δ::Real, accrate::Real)

Update of the random walk step size for the Metropolis-Hastings algorithm
"""
function update_stepsize(δ::Real, accrate::Real)
    Δδ = 0.01 * (2 * (accrate > 0.44) - 1)
    return δ * exp(Δδ)
end