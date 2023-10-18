using PMP

using Distributions, JLD2, SpecialFunctions, Test, CSV, DataFrames, Dates, RollingFunctions

@testset "ExtendedExtremes.jl" begin
    include("distributions/pearsontype1.jl")
    include("distributions/pearsontype1b.jl")
    include("moisture_maximization.jl")
    include("other_PMP_methods.jl")
    include("data_example.jl")
end;