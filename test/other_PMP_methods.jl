
@testset "Other PMP methods" begin
    
    rain = load("data/mm_rain_data.jld2", "Rain")
    date = load("data/mm_rain_data.jld2", "Date")

    # GEV
    @testset "PMP_GEV" begin
        @test PMP.PMP_GEV(rain, date, 100000, 24, 24) == 200.60269271773075
        @test PMP.PMP_GEV(rain, date, 100000) == 200.60269271773075
    end

    # Hershfield
    @testset "PMP_Hershfield" begin      
        @test PMP.PMP_Hershfield(rain, date, 15, 24, 24) ≈ 265.85822091436887  
        @test PMP.PMP_Hershfield(rain, date, 15) ≈ 265.85822091436887
        @test PMP.PMP_Hershfield(rain, date, 24, 24)[1] == 84.41003515317207
        @test PMP.PMP_Hershfield(rain, date)[1] == 84.41003515317207
    end

end

