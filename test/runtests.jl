using PMP

using Distributions, JLD2, SpecialFunctions, Test

@testset "ExtendedExtremes.jl" begin
    include("distributions/pearsontype1.jl")
end;