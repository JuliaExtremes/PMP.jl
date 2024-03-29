@testset "total precipitation" begin

    @testset "maximum rain on d₂ when d₁ ≥ 24h" begin
        rain1 = load("data/mm_shortdata1.jld2", "Rain")
        date1 = load("data/mm_shortdata1.jld2", "Date")

        storm_a = total_precipitation(rain1, date1, 24, 72)
        storm_b = total_precipitation(rain1, DateTime.(date1), 24, 72)

        @test storm_a == storm_b
        @test length(storm_a.Rain) == 4
        @test maximum(storm_a.Rain) == 25.5
    end

    @testset "maximum rain on d₂ when d₁ < 24h" begin
        rain2 = load("data/mm_shortdata2.jld2", "Rain")
        date2 = load("data/mm_shortdata2.jld2", "Date")

        storm_c = total_precipitation(rain2, date2, 3, 24)

        @test storm_c.Date[1] == DateTime(2011, 06, 01, 03)
        @test minimum(storm_c.Rain) == 3.4
    end

end



@testset "storm selection" begin

    rain = load("data/mm_rain_data.jld2", "Rain")
    date = load("data/mm_rain_data.jld2", "Date")

    storm_a = storm_selection(rain, date, 0.1, 24, 72)
    storm_b = storm_selection(rain, date, 0.1)

    @testset "1953-2012" begin
        storm_c = storm_selection(rain, date, 0.1, 24)
        storm_d = storm_selection(rain, date, 0.1)

        @test storm_a == storm_c
        @test storm_b == storm_d
        @test maximum(storm_a.Rain) == 99.5
    end

    @testset "1996" begin
        data1996 = filter(r -> Dates.Year.(r.date) == Year(1996), DataFrames.DataFrame(date = date, rain = rain))

        storm_e = filter(r -> Dates.Year.(r.Date) == Year(1996), storm_a)
        precip1996 = total_precipitation(data1996.rain, data1996.date, 24, 72)

        @test length(storm_e.Rain) == floor(Int, length(precip1996.Rain)/10)+1
        @test maximum(storm_e.Rain) == maximum(precip1996.Rain)
    end

end



@testset "storm PW" begin
    dew = load("data/mm_dew_data.jld2", "Dew")
    date = load("data/mm_dew_data.jld2", "Date")

    @testset "get_maximum_persisting_dew" begin
        @test get_max_persisting_dew(dew[1:48], 1, 12) == 6.4
        @test get_max_persisting_dew(dew[1:48], 3, 24) == 9.0
    end

    @testset "dewpoint_to_PW" begin
        @test dewpoint_to_PW(-10) == 8
        @test dewpoint_to_PW(33) == 123
        @test dewpoint_to_PW(6.4) == 15.4
    end

    @testset "PW_storm" begin
        @test PW_storm([date[1]], dew, DateTime.(date), 48, 1, 12).Dew[] == 6.4
        @test PW_storm([date[1]], dew, date, 48, 1, 12).PW[] == 15.4
        @test Date(PW_storm(DateTime.([date[1]]), dew, date, 48, 1, 12).Date[]) == date[1]
    end
end



@testset "Maximum precipitable water" begin

    pw = load("data/mm_pw_data.jld2", "PW")
    date = load("data/mm_pw_data.jld2", "Date")
    
    @testset "PW_max" begin
        pw_max = PW_max(pw, date)

        @test pw_max[4,2] == 93.6
    end

    @testset "PW_return_period" begin
        pw_rp = PW_return_period(pw, date, 100)

        @test pw_rp[4,2] == 93.74137033150302
    end
end



@testset "Storm maximization" begin
    rain = load("data/mm_rain_data.jld2", "Rain")[10660:end]
    date = load("data/mm_rain_data.jld2", "Date")[10660:end]
    pw = load("data/mm_pw_data.jld2", "PW")[10673:end]
    pw_max = PW_max(pw, date).PW_max

    @testset "Storm" begin
        storm = storm_maximization(rain, pw, date, pw_max)

        @test maximum(storm.Maximized_Rain) == 130.07567567567565 
    end

    @testset "PMP" begin
        pmp = PMP_mm(rain, pw, date, pw_max)
        
        @test pmp == 130.07567567567565
    end
end