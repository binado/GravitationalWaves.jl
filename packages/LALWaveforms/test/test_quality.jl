using TestItems

@testsnippet Setup_quality begin
    using AbstractWaveformsTests

    using Aqua
    using JET
end

@testitem "Aqua" setup = [Common, Setup_quality] begin
    Aqua.test_all(LALWaveforms)
end

@testitem "JET" setup = [Common, Setup_quality] begin
    if JET.JET_AVAILABLE
        @testset "JET" begin
            JET.test_package(LALWaveforms; target_modules = (LALWaveforms,))
        end
    end
end
