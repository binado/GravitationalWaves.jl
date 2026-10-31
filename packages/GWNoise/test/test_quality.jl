using TestItems

@testsnippet Setup_quality begin
    using Aqua
    using JET
end

@testitem "Aqua" setup = [Common, Setup_quality] begin
    Aqua.test_all(GWNoise)
end

@testitem "JET" setup = [Common, Setup_quality] begin
    # JET's analysis is tied to the Julia version it was built for.
    if VERSION.major == 1 && VERSION.minor in (12, 13) && isempty(VERSION.prerelease)
        @testset "JET" begin
            JET.test_package(GWNoise; target_modules = (GWNoise,))
        end
    end
end
