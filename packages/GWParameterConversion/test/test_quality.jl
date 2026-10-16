using Aqua
using JET

@testset "Aqua" begin
    Aqua.test_all(GWParameterConversion)
end

# JET's analysis is tied to the Julia version it was built for.
if VERSION.major == 1 && VERSION.minor == 12 && isempty(VERSION.prerelease)
    @testset "JET" begin
        JET.test_package(GWParameterConversion; target_modules = (GWParameterConversion,))
    end
end
