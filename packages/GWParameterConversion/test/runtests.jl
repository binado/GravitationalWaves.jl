using Test
using GWParameterConversion

@testset "GWParameterConversion" begin
    include("test_masses.jl")
    include("test_spins.jl")
    include("test_convert.jl")
    include("test_mass_coordinates.jl")
    include("test_reference.jl")
    include("test_quality.jl")
    include("test_reactant.jl")
end
