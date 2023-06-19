# modifier les tests : donnees synthetique
@testset "maximum rain on d₂ when d₁ < 24h" begin
    rain = load("data/mm_pwh_data.jld2", "PW")
    date = load("data/mm_pwh_data.jld2", "Date")

    storm = PMP.max_rain(rain, date, 72)
    
    @test storm.Rain[106] == 21.0
    @test maximum(storm.Rain) == 81.0
end



# modifier les tests : donnees synthetique
@testset "maximum rain on d₂ when d₁ ≥ 24h" begin
    rain = load("data/mm_rain_data.jld2", "Rain")
    date = load("data/mm_rain_data.jld2", "Date")

    storm = PMP.max_rain(rain, date, 3, false)
    
    @test storm.Rain[2791] == 15.8
    @test maximum(storm.Rain) == 81.9
end



# modifier les tests : donnees synthetique
@testset "total precipitation" begin
    rain = load("data/mm_rain_data.jld2", "Rain")
    date = load("data/mm_rain_data.jld2", "Date")

    storm = PMP.total_precipitation(rain, date, 24, 72)

    @test storm.Rain[555] == 53.099999999999994

end



# modifier les tests : donnees synthetique
@testset "storm selection" begin
    rain = load("data/mm_rain_data.jld2", "Rain")
    date = load("data/mm_rain_data.jld2", "Date")

    threshold = quantile(rain, 0.95)
    storm = PMP.storm_selection(rain, date, 0.1, 24, 72)

    @test storm[7, :].Rain == 30.700000000000003
    @test storm[7, :].Date == Dates.Date(1954, 06, 25)
    
end



@testset "get_max_persisting_dew" begin
    dew = load("data/mm_dew_data.jld2", "Dew")

    @test PMP.get_max_persisting_dew(dew, 12) == 6.4
end



@testset "dewpoint_to_PW" begin
    @test PMP.dewpoint_to_PW(-10) == 8
    @test PMP.dewpoint_to_PW(33) == 123
    @test PMP.dewpoint_to_PW(10.5) == 22
end



@testset "Maximum precipitable water" begin
    
    pw = load("data/mm_pw_data.jld2", "PW")
    date = load("data/mm_pw_data.jld2", "Date")
    
    @testset "PW_max" begin
        pw_max = PMP.PW_max(pw, date)

        @test pw_max[4,2] == 93.6
    end

    @testset "PW_return_period" begin
        pw_rp = PMP.PW_return_period(pw, date, 100)

        @test pw_rp[4,2] == 93.74137033150302
    end
end



@testset "PMP_mm" begin
    rain = load("data/mm_rain_data.jld2", "Rain")[10660:end]
    date = load("data/mm_rain_data.jld2", "Date")[10660:end]
    pw = load("data/mm_pw_data.jld2", "PW")[10673:end]
    pw_max = PMP.PW_max(pw, date).PW_max

    PMPa, PMPb, storms = PMP.PMP_mm(rain, pw, date, pw_max)

    @test PMPa == 130.07567567567565
    @test PMPb == 112.8
end