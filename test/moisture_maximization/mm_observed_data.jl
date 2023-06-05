@testset "storm_selection_cluster" begin
    rain = load("data/mm_rain_data.jld2", "Rain")
    date = load("data/mm_rain_data.jld2", "Date")
    threshold = quantile(rain, 0.95)
    storm_PMP = PMP.storm_selection_cluster(rain, date, threshold)

    @test storm_PMP[7, :].Duration == 2
    @test storm_PMP[7, :].Rain_max == 39.8
    @test storm_PMP[7, :].Rain_total == 43.0
    @test storm_PMP[7, :].Date == Dates.Date(2011, 08, 21)
end

@testset "storm_selection_fixed" begin
    rain = load("data/mm_rain_data.jld2", "Rain")
    date = load("data/mm_rain_data.jld2", "Date")
    threshold = quantile(rain, 0.95)
    storm_PMP = PMP.storm_selection_fixed(rain, date, threshold)

    @test storm_PMP[7, :].Rain == 39.8
    @test storm_PMP[7, :].Date == Dates.Date(2011, 08, 21)
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



@testset "PW_max" begin
    pw = load("data/mm_pw_data.jld2", "PW")
    date = load("data/mm_pw_data.jld2", "Date")
    pw_max = PMP.PW_max(pw, date)

    @test PMP.PW_max(pw, date)[3,2] == 81.0
end
