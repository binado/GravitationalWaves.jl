using Test
using Aqua
using JET
using AbstractGWWaveforms
using AbstractGWWaveforms: check_frequency_grid, frequency_bounds

include("contract.jl")

struct Dummy <: AbstractWaveformGenerator
    maximum::Union{Nothing, Float64}
end
AbstractGWWaveforms.approximant(::Dummy) = "Dummy"
AbstractGWWaveforms.reference_frequency(::Dummy) = 50.0
AbstractGWWaveforms.minimum_frequency(::Dummy) = 20.0
AbstractGWWaveforms.maximum_frequency(g::Dummy) = g.maximum

function AbstractGWWaveforms.hphc_fd(f::AbstractRange, params::NamedTuple, g::Dummy)
    check_frequency_grid(f)
    foreach(k -> k in (:mass_1, :mass_2) || throw(ArgumentError("unknown parameter $k")), keys(params))
    f_min, f_max = frequency_bounds(g, f)
    h = [f_min <= fi <= f_max ? complex(1.0, fi) : 0.0im for fi in f]
    return (; plus = h, cross = h)
end

@testset "AbstractGWWaveforms" begin
    @testset "Aqua" begin
        Aqua.test_all(AbstractGWWaveforms)
    end

    @testset "JET" begin
        JET.test_package(AbstractGWWaveforms; target_modules = (AbstractGWWaveforms,))
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

    test_generator_contract(Dummy(nothing); params = (mass_1 = 30.0, mass_2 = 30.0), f = 0:0.5:100)
    test_generator_contract(Dummy(60.0); params = (mass_1 = 30.0, mass_2 = 30.0), f = 0:0.5:100)
end
