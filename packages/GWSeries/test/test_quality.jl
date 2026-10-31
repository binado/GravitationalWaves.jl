using TestItems

@testsnippet Setup_quality begin
    using Aqua
    using JET
end

@testitem "Aqua" setup = [Common, Setup_quality] begin
    Aqua.test_all(GWSeries)
end

@testitem "JET" setup = [Common, Setup_quality] begin
    if VERSION.major == 1 && VERSION.minor in (12, 13) && isempty(VERSION.prerelease)
        @testset "JET" begin
            JET.test_package(GWSeries; target_modules = (GWSeries,))
        end
    end
end
