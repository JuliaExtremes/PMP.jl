
@testset "PearsonType1 constructor" begin
    @test PearsonType1(0., 1., 1. ,1.) == PearsonType1{Float64}(0.,1.,1.,1.)
    @test PearsonType1(0, 1, 1, 1) == PearsonType1{Float64}(0.,1.,1.,1.)
    @test PearsonType1(Float32(0), 1, 1, 1) == PearsonType1{Float32}(0.0f0, 1.0f0, 1.0f0, 1.0f0)
    @test  PearsonType1() == PearsonType1{Float64}(0.,1.,1.,1.)
end

@testset "PearsonType1 parameters" begin
    pd = PearsonType1(-1.,1.,2.,3.)
    @test location(pd) ≈ -1.
    @test scale(pd) ≈ 2.
    @test all(shape(pd) .≈ (2., 3.))
    @test all(params(pd) .≈ (-1., 1., 2., 3.))
end

@testset "PearsonType1 evaluations" begin
    pd = PearsonType1(-1.,1.,2.,3.)
    x = -.5
    
    @testset "cdf" begin
        @test cdf(pd, x) ≈ cdf(Beta(shape(pd)...),(x-location(pd))/scale(pd))
    end
    
    @testset "insupport" begin
        @test insupport(pd, x)
        @test !insupport(pd, 2.)
    end
    
    @testset "logpdf" begin
        a, b, α, β = params(pd)
        true_log_pdf_at_x = -SpecialFunctions.logbeta(α, β) + (α-1)*log(x-a) + (β-1)*log(b-x) - (α+β-1)*log(b-a)
        
        @test logpdf(pd, x) ≈ true_log_pdf_at_x
    end
    
    @testset "maximum" begin
       @test maximum(pd) == params(pd)[2] 
    end
    
    @testset "minimum" begin
       @test minimum(pd) == params(pd)[1] 
    end
    
    @testset "quantile" begin
        @test quantile(pd, .9) ≈ location(pd) + scale(pd) * quantile(Beta(shape(pd)...), .9)
    end
    
    @testset "rand" begin
        @test rand(pd) isa Any
    end
end



@testset "PearsonType1 statistics" begin
    pd = PearsonType1(-1.,1.,2.,3.)
    α = 2.
    β = 3.
    
    x = -.5
    
    @testset "mean" begin
        @test mean(pd) == scale(pd)*α/(α+β) + location(pd)
    end

    @testset "var" begin
        @test var(pd) == scale(pd)^2*α*β/((α+β)^2*(α+β+1))
    end
    
    @testset "std" begin
        @test std(pd) == sqrt(var(pd))
    end

    @testset "modes" begin
        @test modes(pd) == location(pd) .+ scale(pd) .* modes(Beta(shape(pd)...))
    end
    
    @testset "mode" begin
        @test mode(pd) == location(pd) + scale(pd)* mode(Beta(shape(pd)...))
    end

    @testset "skewness" begin
        @test skewness(pd) == 2*(β-α)*sqrt(α+β+1)/((α+β+2)*sqrt(α*β))
    end
    
    @testset "kurtosis" begin
        @test kurtosis(pd) == 6*(α^3-α^2*(2*β-1)+β^2*(β+1)-2*α*β*(β+2))/(α*β*(α+β+2)*(α+β+3))
        @test kurtosis(pd, false) ≈ kurtosis(Beta(shape(pd)...)) + 3
        @test kurtosis(pd, false) == kurtosis(pd) + 3
        @test kurtosis(pd, true) == kurtosis(pd)                                     
    end

    @testset "entropy" begin
        base = 3
        @test entropy(pd) ≈ entropy(Beta(shape(pd)...)) + log(scale(pd))
        @test entropy(pd, base) ≈ entropy(pd)/log(base)
    end
end



@testset "fit_mme PearsonType1" begin
    y = load("data/pearsontype1_sample.jld2", "y")
    a = -1
    fd1 = fit_mme(PearsonType1, y, a)
    b = fd1.b
    α, β = shape(fd1)

    @test mean(y) ≈ (b-a)*α/(α+β)+a atol=0.05
    @test var(y) ≈ (b-a)^2*α*β/((α+β)^2*(α+β+1)) atol=0.01
    @test skewness(y) ≈ 2*(β-α)*sqrt(α+β+1)/((α+β+2)*sqrt(α*β))
    @test kurtosis(y) ≈ 6*(α^3-α^2*(2*β-1)+β^2*(β+1)-2*α*β*(β+2))/(α*β*(α+β+2)*(α+β+3))

    fd2 = fit_mme(PearsonType1, y)
    a, b, α, β = params(fd2)

    @test mean(y) ≈ (b-a)*α/(α+β)+a atol=0.01
    @test var(y) ≈ (b-a)^2*α*β/((α+β)^2*(α+β+1)) atol=0.01
    @test skewness(y) ≈ 2*(β-α)*sqrt(α+β+1)/((α+β+2)*sqrt(α*β))
    @test kurtosis(y) ≈ 6*(α^3-α^2*(2*β-1)+β^2*(β+1)-2*α*β*(β+2))/(α*β*(α+β+2)*(α+β+3))
end



@testset "fit_mle PearsonType1" begin
    y = load("data/pearsontype1_sample.jld2", "y")
    a = -1
    fd1 = fit_mle(PearsonType1, y, [maximum(y), 1., 2.], a)

    @test maximum(fd1) ≈ 1. atol=.05
    @test shape(fd1)[1] ≈ 2. atol=.1
    @test shape(fd1)[2] ≈ 3. atol=.3

    fd2 = fit_mle(PearsonType1, y, [minimum(y), maximum(y), 1., 2.])
    
    @test minimum(fd2) ≈ -1. atol=.1
    @test maximum(fd2) ≈ 1. atol=.1
    @test shape(fd2)[1] ≈ 2. atol=.1
    @test shape(fd2)[2] ≈ 3. atol=.3

    fd3 = fit_mle(PearsonType1, y, a)

    @test maximum(fd3) ≈ 1. atol=.05
    @test shape(fd3)[1] ≈ 2. atol=.1
    @test shape(fd3)[2] ≈ 3. atol=.3

    fd4 = fit_mle(PearsonType1, y)

    @test minimum(fd4) ≈ -1. atol=.01
    @test maximum(fd4) ≈ 1. atol=.1
    @test shape(fd4)[1] ≈ 2. atol=.3
    @test shape(fd4)[2] ≈ 3. atol=.3
end



@testset "getinitialvalues PearsonType1" begin
    y = load("data/pearsontype1_sample.jld2", "y")
    a = -1
    ivalues1 = getinitialvalues(PearsonType1, y, a)

    @test ivalues1[1] == -1.
    @test ivalues1[2] ≈ 1. atol=.1
    @test ivalues1[3] ≈ 2. atol=.3
    @test ivalues1[4] ≈ 3. atol=.4

    ivalues2 = getinitialvalues(PearsonType1, y)

    @test ivalues2[1] ≈ -1. atol=.01
    @test ivalues2[2] ≈ 1. atol=.1
    @test ivalues2[3] ≈ 2. atol=.3
    @test ivalues2[4] ≈ 3. atol=.4
end