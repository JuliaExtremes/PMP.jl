@testset "data.jl" begin
    
    @testset "dataset(name)" begin
        # nonexistent file throws
        @test_throws ErrorException PMP.dataset("nonexistant")

        # rain loading
        @test_logs df = PMP.dataset("rain")
    end

end