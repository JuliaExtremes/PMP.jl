
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
    dist = PMP.getdistribution(pd)
    
    x = -.5
    
    @testset "cdf" begin
        @test cdf(pd, x) ≈ cdf(dist, x)
    end
    
    @testset "getdistribution" begin

        p = params(dist)

        @test p[1] ≈ location(pd)
        @test p[2] ≈ scale(pd)
        @test p[3] == Beta(shape(pd)...)

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
        @test quantile(pd, .9) ≈ quantile(dist, .9)
    end
    
    @testset "rand" begin
        @test rand(pd) isa Any
    end
end



@testset "PearsonType1 statistics" begin
    pd = PearsonType1(-1.,1.,2.,3.)
    dist = PMP.getdistribution(pd)
    α = 2.
    β = 3.
    
    x = -.5
    
    @testset "mean" begin
        @test mean(pd) ≈ mean(dist)
        @test mean(pd) == scale(pd)*α/(α+β) + location(pd)
    end

    @testset "var" begin
        @test var(pd) ≈ var(dist)
        @test var(pd) == scale(pd)^2*α*β/((α+β)^2*(α+β+1))
    end
    
    @testset "std" begin
        @test std(pd) ≈ std(dist)
        @test std(pd) == sqrt(var(pd))
    end

    @testset "modes" begin
        @test modes(pd) ≈ modes(dist)
    end
    
    @testset "mode" begin
        @test mode(pd) ≈ mode(dist)
    end

    @testset "skewness" begin
        @test skewness(pd) ≈ skewness(dist)
        @test skewness(pd) == 2*(β-α)*sqrt(α+β+1)/((α+β+2)*sqrt(α*β))
    end
    
    @testset "kurtosis" begin
        @test kurtosis(pd) ≈ kurtosis(dist)
        @test kurtosis(pd) == 6*(α^3-α^2*(2*β-1)+β^2*(β+1)-2*α*β*(β+2))/(α*β*(α+β+2)*(α+β+3))
        @test kurtosis(pd, false) ≈ kurtosis(dist, false)
        @test kurtosis(pd, false) == kurtosis(pd) + 3                                     
    end

    @testset "entropy" begin
        base = 3
        @test entropy(pd) ≈ entropy(dist)
        @test entropy(pd, base) ≈ entropy(dist, base)
    end
end



@testset "fit_mme" begin
    y = load("data/persontype1_sample.jld2", "y")
    fd = PMP.fit_mme(PearsonType1, y)
    a, b, α, β = params(fd)

    @test mean(y) ≈ (b-a)*α/(α+β) + a
    @test var(y) ≈ (b-a)^2*α*β/((α+β)^2*(α+β+1))
    @test skewness(y) ≈ 2*(β-α)*sqrt(α+β+1)/((α+β+2)*sqrt(α*β))
    @test kurtosis(y) ≈ 6*(α^3-α^2*(2*β-1)+β^2*(β+1)-2*α*β*(β+2))/(α*β*(α+β+2)*(α+β+3))
end



#@testset "fit_mle" begin
#    y = load("data/persontype1_sample.jld2", "y")
#    fd = PMP.fit_mle(PearsonType1, y, [minimum(y), maximum(y), 1., 1.])
    
 #   @test minimum(fd) ≈ -1. atol=0.01
#    @test maximum(fd) ≈ 1. atol=0.05
#    @test shape(fd)[1] ≈ 2. atol=.1
#    @test shape(fd)[2] ≈ 3. atol=.3

    #fd2 = PMP.fit_mle(PearsonType1, y)

    #@test minimum(fd) ≈ -1. atol=0.01
    #@test maximum(fd) ≈ 1. atol=0.05
    #@test shape(fd)[1] ≈ 2. atol=.1
    #@test shape(fd)[2] ≈ 3. atol=.3
#end