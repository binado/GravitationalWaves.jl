using Test
using GWWaveforms
using Enzyme: Enzyme
using Mooncake: Mooncake

include("../../packages/GWWaveforms/test/ad_helpers.jl")

test_waveform_gradients(
    Pair{String, Any}[
        "ForwardDiff" => AutoForwardDiff(),
        "Mooncake" => AutoMooncake(; config = nothing),
        "Enzyme" => AutoEnzyme(; function_annotation = Enzyme.Const),
    ]
)

include("test_reactant.jl")
