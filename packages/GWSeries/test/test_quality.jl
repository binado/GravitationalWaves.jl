using Aqua
using JET

@testset "Aqua" begin
    Aqua.test_all(GWSeries)
end

if VERSION.major == 1 && VERSION.minor == 12 && isempty(VERSION.prerelease)
    @testset "JET" begin
        JET.test_package(GWSeries; target_modules = (GWSeries,))
    end
end
