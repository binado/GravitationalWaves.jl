include("ad_helpers.jl")

test_waveform_gradients(["ForwardDiff" => AutoForwardDiff()])
