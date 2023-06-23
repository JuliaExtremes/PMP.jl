
@testset "Other PMP methods" begin
    
    rain = load("data/mm_rain_data.jld2", "Rain")
    date = load("data/mm_rain_data.jld2", "Date")

    # GEV
    @testset "PMP_GEV" begin
        @test PMP.PMP_GEV(rain, date, 100000) == 200.60269271773075
    end

    # Hershfield
    @testset "PMP_Hershfield" begin        
        @test PMP.PMP_Hershfield(rain, date, 15) == 265.8582209143689
        @test PMP.PMP_Hershfield(rain, date)[1] == 84.41003515317207
    end

end

