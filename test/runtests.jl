using PMP

using CSV, DataFrames, Distributions, SpecialFunctions, Test

@testset "ExtendedExtremes.jl" begin
    include("distributions/pearsontype1.jl")
end;