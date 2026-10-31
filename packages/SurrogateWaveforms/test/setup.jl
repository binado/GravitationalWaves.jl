using TestItems

# Evaluated afresh in each test item, including mutable fixtures and dummy types.
@testsnippet Common begin
    using GWSeries
    using Test
    using SurrogateWaveforms

    const FIXTURE_DIR = joinpath(@__DIR__, "..", "..", "..", "reference", "surrogates")
    const MODEL_PARAMS = (NRSur7dq4 = SurrogateWaveforms.NRSur7dq4Params(), NRSur7dq4v2 = SurrogateWaveforms.NRSur7dq4v2Params())
end
