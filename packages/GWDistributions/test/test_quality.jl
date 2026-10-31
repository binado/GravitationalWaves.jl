using TestItems

@testsnippet Setup_quality begin
    using Aqua
    using JET
end

@testitem "Aqua" setup = [Common, Setup_quality] begin
    Aqua.test_all(GWDistributions)
end

@testitem "JET" setup = [Common, Setup_quality] begin
    # JET's analysis is tied to the Julia version it was built for.
    if VERSION.major == 1 && VERSION.minor == 12 && isempty(VERSION.prerelease)
        @testset "JET" begin
            JET.test_package(GWDistributions; target_modules = (GWDistributions,))
        end
    end
end
