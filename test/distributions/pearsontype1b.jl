
@testset "PearsonType1b constructor" begin
    @test PearsonType1b(1., 1. ,1.) == PearsonType1b{Float64}(1.,1.,1.)
    @test PearsonType1b(1, 1, 1) == PearsonType1b{Float64}(1.,1.,1.)
    @test PearsonType1b(1, 1, 1) == PearsonType1b{Float32}(1.0f0, 1.0f0, 1.0f0)
    @test  PearsonType1b() == PearsonType1b{Float64}(1.,1.,1.)
end

@testset "PearsonType1b parameters" begin
    pd = PearsonType1b(1.,2.,3.)
    @test location(pd) == 0
    @test scale(pd) ≈ 1.
    @test all(shape(pd) .≈ (2., 3.))
    @test all(params(pd) .≈ (1., 2., 3.))
end

@testset "PearsonType1b evaluations" begin
    pd = PearsonType1b(1.,2.,3.)
    x = .5
    
    @testset "cdf" begin
        @test cdf(pd, x) ≈ cdf(Beta(shape(pd)...), x/pd.b)
    end
    
    @testset "insupport" begin
        @test insupport(pd, x)
        @test !insupport(pd, 2.)
    end
    
    @testset "logpdf" begin
        b, α, β = params(pd)
        true_log_pdf_at_x = -SpecialFunctions.logbeta(α, β) + (α-1)*log(x) + (β-1)*log(b-x) - (α+β-1)*log(b)
        
        @test logpdf(pd, x) ≈ true_log_pdf_at_x
    end
    
    @testset "maximum" begin
       @test maximum(pd) == scale(pd)
    end
    
    @testset "minimum" begin
       @test minimum(pd) == 0
    end
    
    @testset "quantile" begin
        @test quantile(pd, .9) ≈ pd.b * quantile(Beta(shape(pd)...), .9)
    end
    
    @testset "rand" begin
        @test rand(pd) isa Any
    end
end



@testset "PearsonType1b statistics" begin
    pd = PearsonType1b(1.,2.,3.)
    α = 2.
    β = 3.
    
    x = .5
    
    @testset "mean" begin
        @test mean(pd) == scale(pd)*α/(α+β)
    end

    @testset "var" begin
        @test var(pd) == scale(pd)^2*α*β/((α+β)^2*(α+β+1))
    end
    
    @testset "std" begin
        @test std(pd) == sqrt(var(pd))
    end

    @testset "modes" begin
        @test modes(pd) == scale(pd) .* modes(Beta(shape(pd)...))
    end
    
    @testset "mode" begin
        @test mode(pd) == scale(pd)* mode(Beta(shape(pd)...))
    end

    @testset "skewness" begin
        @test skewness(pd) == 2*(β-α)*sqrt(α+β+1)/((α+β+2)*sqrt(α*β))
    end
    
    @testset "kurtosis" begin
        @test kurtosis(pd) == 6*(α^3-α^2*(2*β-1)+β^2*(β+1)-2*α*β*(β+2))/(α*β*(α+β+2)*(α+β+3))
        @test kurtosis(pd, false) ≈ kurtosis(Beta(shape(pd)...)) + 3
        @test kurtosis(pd, false) == kurtosis(pd) + 3                                     
    end

    @testset "entropy" begin
        base = 3
        @test entropy(pd) ≈ entropy(Beta(shape(pd)...)) + log(scale(pd))
        @test entropy(pd, base) ≈ entropy(pd)/log(base)
    end
end



@testset "fit_mme" begin
    y = load("data/pearsontype1b_sample.jld2", "y")
    fd = PMP.fit_mme(PearsonType1b, y)
    b, α, β = params(fd)

    @test var(y) ≈ b^2*α*β/((α+β)^2*(α+β+1))
    @test skewness(y) ≈ 2*(β-α)*sqrt(α+β+1)/((α+β+2)*sqrt(α*β))
    @test kurtosis(y) ≈ 6*(α^3-α^2*(2*β-1)+β^2*(β+1)-2*α*β*(β+2))/(α*β*(α+β+2)*(α+β+3))
end