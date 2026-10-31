using TestItems

@testsnippet Setup_ad begin
    using AbstractWaveformsTests

    include("ad_helpers.jl")
end

@testitem "Waveform parameter gradients" tags = [:ad] setup = [Common, Setup_ad] begin
    test_waveform_gradients(["ForwardDiff" => AutoForwardDiff()])
end
