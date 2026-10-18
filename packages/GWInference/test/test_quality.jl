using Aqua
using JET

@testset "Aqua" begin
    Aqua.test_all(GWInference)
end

if VERSION.major == 1 && VERSION.minor == 12 && isempty(VERSION.prerelease)
    @testset "JET" begin
        JET.test_package(GWInference; target_modules = (GWInference,))
    end
end
