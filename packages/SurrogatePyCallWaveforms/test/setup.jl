using TestItems

# Evaluated afresh in each test item, including mutable fixtures and dummy types.
@testsnippet Common begin
    using Test
    using PyCall
    import AbstractWaveforms
    import SurrogateWaveforms
    using AbstractWaveformsTests: test_modes_td_contract
    using SurrogatePyCallWaveforms

    # Needs a Python with gwsurrogate (PyCall built against it) and the NRSur7dq4 model data. Without
    # them the tests are skipped locally and fail on CI.
    function python_ready()
        return try
            pyimport("gwsurrogate.surrogate")
            true
        catch err
            err isa InterruptException && rethrow()
            false
        end
    end
end
