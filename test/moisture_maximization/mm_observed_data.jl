@testset "storm selection" begin

    rain = load("data/mm_rain_data.jld2", "Rain")
    date = load("data/mm_rain_data.jld2", "Date")
    threshold = quantile(rain, 0.95)

    @testset "storm_selection_cluster" begin
        storm_PMP = PMP.storm_selection_cluster(rain, date, threshold)

        @test storm_PMP[7, :].Duration == 2
        @test storm_PMP[7, :].Rain_max == 39.8
        @test storm_PMP[7, :].Rain_total == 43.0
        @test storm_PMP[7, :].Date == Dates.Date(2011, 08, 21)
    end

    @testset "storm_selection_fixed" begin
        storm_PMP = PMP.storm_selection_fixed(rain, date, threshold)

        @test storm_PMP[7, :].Rain == 39.8
        @test storm_PMP[7, :].Date == Dates.Date(2011, 08, 21)
    end
end



@testset "get_max_persisting_dew" begin
    dew = load("data/mm_dew_data.jld2", "Dew")

    @test PMP.get_max_persisting_dew(12, dew) == 6.4
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
        pw_rp = PMP.PW_return_period(100, pw, date)

        @test pw_rp[4,2] == 93.74137033150302
    end
end



@testset "PMP_mm" begin
    rain = load("data/mm_rain_data.jld2", "Rain")
    date = load("data/mm_rain_data.jld2", "Date")
    pw = load("data/mm_pw_data.jld2", "PW")[10673:end]
    pw_max = PMP.PW_max(pw, date).PW_max

    PMPa, PMPb, storms = PMP.PMP_mm(rain, pw, date, pw_max)

    @test PMPa == 130.07567567567565
    @test PMPb == 112.8
end