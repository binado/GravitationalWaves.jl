using GWSeries
using Test
using Aqua
using JET
using AbstractWaveforms
using AbstractWaveforms: check_frequency_grid, frequency_bounds
using AbstractWaveforms: AlignedSpinQuadrupolarWaveformGenerator, amplitude_phase_fd
using AbstractWaveformsTests

struct Dummy <: AbstractWaveformGenerator
    maximum::Union{Nothing, Float64}
end
AbstractWaveforms.approximant(::Dummy) = "Dummy"
AbstractWaveforms.reference_frequency(::Dummy) = 50.0
AbstractWaveforms.minimum_frequency(::Dummy) = 20.0
AbstractWaveforms.maximum_frequency(g::Dummy) = g.maximum

function AbstractWaveforms.hphc_fd(f::AbstractRange{Float64}, params::NamedTuple, g::Dummy)
    check_frequency_grid(f)
    f_min, f_max = frequency_bounds(g, f)
    h = [f_min <= fi <= f_max ? complex(1.0, fi) : 0.0im for fi in f]
    return hcat(h, h)
end

function AbstractWaveforms.validate_parameters(f::AbstractRange{Float64}, params::NamedTuple, ::Dummy)
    foreach(k -> k in (:mass_1, :mass_2) || throw(ArgumentError("unknown parameter $k")), keys(params))
    return nothing
end

# Implements only amplitude_phase_fd; polarization assembly comes from the interface.
struct DummyQuadrupole <: AlignedSpinQuadrupolarWaveformGenerator end
AbstractWaveforms.approximant(::DummyQuadrupole) = "DummyQuadrupole"
AbstractWaveforms.reference_frequency(::DummyQuadrupole) = 50.0
AbstractWaveforms.minimum_frequency(::DummyQuadrupole) = 20.0
AbstractWaveforms.maximum_frequency(::DummyQuadrupole) = 60.0

function AbstractWaveforms.amplitude_phase_fd(f::AbstractRange{Float64}, params::NamedTuple, g::DummyQuadrupole)
    check_frequency_grid(f)
    lo, hi = frequency_bounds(g, f)
    amplitude = [lo <= fi <= hi ? 2.0 : 0.0 for fi in f]
    phase = [lo <= fi <= hi ? 0.1 * fi - 2get(params, :phase, 0.0) : 0.0 for fi in f]
    return (; amplitude, phase)
end

function AbstractWaveforms.validate_parameters(f::AbstractRange{Float64}, params::NamedTuple, ::DummyQuadrupole)
    check_frequency_grid(f)
    foreach(k -> k == :theta_jn || throw(ArgumentError("unknown parameter $k")), keys(params))
    return nothing
end

include("test_modes.jl")

@testset "AbstractWaveforms" begin
    @testset "Aqua" begin
        Aqua.test_all(AbstractWaveforms)
    end

    @testset "JET" begin
        JET.test_package(AbstractWaveforms; target_modules = (AbstractWaveforms,))
    end

    @testset "capability fallbacks" begin
        g = Dummy(nothing)
        @test spin_support(g) === :aligned
        @test supports_mode_array(g) === false
        @test is_differentiable(g) === false
    end

    @testset "check_frequency_grid" begin
        @test check_frequency_grid(0:0.5:10) == 0:0.5:10
        @test_throws ArgumentError check_frequency_grid(1:0.5:10)
        @test_throws ArgumentError check_frequency_grid(0:-0.5:-10)
        @test_throws ArgumentError check_frequency_grid(range(0, 0; length = 1))
    end

    @testset "frequency_bounds" begin
        f = 0:0.5:100
        @test frequency_bounds(Dummy(nothing), f) == (20.0, 100.0)
        @test frequency_bounds(Dummy(60.0), f) == (20.0, 60.0)
    end

    @testset "show" begin
        @test sprint(show, Dummy(60.0)) ==
            "Dummy(\"Dummy\"; reference_frequency=50.0, minimum_frequency=20.0, maximum_frequency=60.0)"
    end

    test_generator_contract(DummyQuadrupole(); params = (; theta_jn = 0.4), f = 0:0.5:100)

    test_generator_contract(Dummy(nothing); params = (mass_1 = 30.0, mass_2 = 30.0), f = 0:0.5:100)
    test_generator_contract(Dummy(60.0); params = (mass_1 = 30.0, mass_2 = 30.0), f = 0:0.5:100)
end
