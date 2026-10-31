using TestItems

@testsnippet Setup_quality begin
    using Aqua
    using JET
end

@testitem "Aqua" setup = [Common, Setup_quality] begin
    Aqua.test_all(GWDistributions)
end

@testitem "JET" setup = [Common, Setup_quality] begin
    if JET.JET_AVAILABLE
        @testset "JET" begin
            JET.test_package(GWDistributions; target_modules = (GWDistributions,))
        end
    end
end
