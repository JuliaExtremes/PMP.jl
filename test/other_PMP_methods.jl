
@testset "Other PMP methods" begin
    
    rain = load("data/mm_rain_data.jld2", "Rain")
    date = load("data/mm_rain_data.jld2", "Date")

    # GEV
    @testset "PMP_GEV" begin
        PMP = PMP.PMP_GEV(rain, date, 100000)

        @test PMP[1] == 200.60269271773075
    end

    # Hershfield
    @testset "PMP_Hershfield" begin
        PMPk = PMP.PMP_Hershfield(rain, date, 15)
        PMPsk = PMP.PMP_Hershfield(rain, date)
        
        @test PMPk == 265.8582209143689
        @test PMPsk[1] == 84.41003515317207
    end
end

