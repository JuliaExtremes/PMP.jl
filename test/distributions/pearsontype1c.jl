
@testset "PearsonType1c constructor" begin
    @test PearsonType1c(10., .5, 1.) == PearsonType1c{Float64}(10.,.5,1.)
    @test PearsonType1c(10, .5, 1) == PearsonType1c{Float64}(10.,.5,1.)
    @test PearsonType1c(10, .5, 1) == PearsonType1c{Float32}(10.0f0, .5f0, 1.0f0)
    @test  PearsonType1c() == PearsonType1c{Float64}(1.,0.5,2.)
end



@testset "PearsonType1c parameters" begin
    pd = PearsonType1c(10., .4, 3.)
    @test location(pd) == 0
    @test scale(pd) ≈ 10.
    @test all(shape(pd) .≈ (.4, 3.))
    @test all(params(pd) .≈ (10., .4, 3.))
end



@testset "PearsonType1c evaluations" begin
    pd = PearsonType1c(10., 0.4, 3.)
    x = .5
    
    @testset "logpdf" begin
        b, μ, ν = params(pd)
        true_log_pdf_at_x = loggamma(ν) - loggamma(μ*ν) - loggamma(ν*(1. - μ)) + (μ*ν - 1.)*log(x) + (ν*(1. - μ) - 1.)*log(b - x) - (ν - 1.)*log(b)
        
        @test logpdf(pd, x) ≈ true_log_pdf_at_x
    end
end



@testset "fit_mle PearsonType1c" begin
    y = load("data/pearsontype1b_sample.jld2", "y")
    fd = fit_mle(PearsonType1c, y, [maximum(y), 1/2, 3.])
    
    @test scale(fd) ≈ 1. atol=.1
    @test shape(fd)[1] ≈ 1/3 atol=.1
    @test shape(fd)[2] ≈ 5. atol=.1

    fd2 = fit_mle(PearsonType1c, y)

    @test scale(fd2) ≈ 1. atol=.1
    @test shape(fd2)[1] ≈ 1/3 atol=.1
    @test shape(fd2)[2] ≈ 5. atol=.1

    x = [Inf, 3]
    fd3 = fit_mle(PearsonType1c, x, [10, 0.5, 3])

    @test_logs (:warn, "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values.")
    @test fd3 == PearsonType1c(10, 0.5, 3)
end



@testset "getinitialvalues PearsonType1c" begin
    y = load("data/pearsontype1b_sample.jld2", "y")
    ivalues = getinitialvalues(PearsonType1c, y)

    @test ivalues[1] ≈ 1. atol=.1
    @test ivalues[2] ≈ 1/3 atol=.1
    @test ivalues[3] ≈ 5. atol=.1
end



@testset "fit_bayes PearsonType1c" begin
    y = load("data/pearsontype1b_sample.jld2", "y")
    trace = fit_bayes(PearsonType1c, y, 1, 1000, 200)

    @test mean(trace[1]) ≈ 1. atol=.01
    @test mean(trace[2]) ≈ 1/3 atol=.1
    @test mean(trace[3]) ≈ 5. atol=.1
end