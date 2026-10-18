using Aqua
using JET

@testset "Aqua" begin
    Aqua.test_all(GWMatchedFiltering)
end

if VERSION.major == 1 && VERSION.minor == 12 && isempty(VERSION.prerelease)
    @testset "JET" begin
        JET.test_package(GWMatchedFiltering; target_modules = (GWMatchedFiltering,))
    end
end
