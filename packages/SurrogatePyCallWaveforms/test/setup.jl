using TestItems

# Evaluated afresh in each test item, including mutable fixtures and dummy types.
@testsnippet Common begin
    using Test
    using PyCall
    import AbstractWaveforms
    import SurrogateWaveforms
    using AbstractWaveformsTests: test_modes_td_contract
    using SurrogatePyCallWaveforms

    # Only checks that PyCall can import gwsurrogate; model data is resolved later. When this is false
    # the `ModelFixture` warns (errors with CI=true) and the `if READY` bodies record `@test_skip`.
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
