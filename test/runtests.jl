using PMP

using Distributions, JLD2, SpecialFunctions, Test, CSV, DataFrames, Dates

@testset "ExtendedExtremes.jl" begin
    include("distributions/pearsontype1.jl")
    include("moisture_maximization/mm_observed_data.jl")
end;