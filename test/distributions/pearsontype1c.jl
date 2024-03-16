
@testset "PearsonType1c constructor" begin
    @test PearsonType1c(10., 1., 1.) == PearsonType1c{Float64}(10.,1.,1.)
    @test PearsonType1c(10, 1, 1) == PearsonType1c{Float64}(10.,1.,1.)
    @test PearsonType1c(10, 1, 1) == PearsonType1c{Float32}(10.0f0, 1.0f0, 1.0f0)
    @test  PearsonType1c() == PearsonType1c{Float64}(1.,0.5,2.)
end

@testset "PearsonType1c parameters" begin
    pd = PearsonType1c(10.,2.,3.)
    @test location(pd) == 0
    @test scale(pd) ≈ 10.
    @test all(shape(pd) .≈ (2., 3.))
    @test all(params(pd) .≈ (10., 2., 3.))
end

@testset "PearsonType1b evaluations" begin
    pd = PearsonType1c(10.,2.,3.)
    x = .5
    
    @testset "logpdf" begin
        b, μ, ν = params(pd)
        true_log_pdf_at_x = loggamma(ν) - loggamma(μ*ν) - loggamma(ν*(1 - μ)) + (μ*ν - 1)*log(x) + (ν*(1 - μ) - 1)*log(b - x) - (ν - 1)*log(b)
        
        @test logpdf(pd, x) ≈ true_log_pdf_at_x
    end
end