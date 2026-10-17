using Test

# Shared checks for every `AbstractWaveformGenerator`. Backends `include` this file from their
# tests, with the interface (`hphc_fd`, the accessors, ...) in scope.
#
# `params` must be valid for `g`, `f` a valid grid, and the support of `g` must lie inside `f`
# with some of it strictly inside, so the masking can be checked on both sides.
function test_generator_contract(g; params, f)
    return @testset "contract: $(sprint(show, g))" begin
        @testset "accessors" begin
            @test approximant(g) isa String
            @test reference_frequency(g) isa Real
            @test minimum_frequency(g) isa Real
            @test maximum_frequency(g) isa Union{Nothing, Real}
            @test spin_support(g) in (:none, :aligned, :precessing)
            @test supports_mode_array(g) isa Bool
            @test is_differentiable(g) isa Bool
            @test occursin(approximant(g), sprint(show, g))
        end

        @testset "grid validation" begin
            @test_throws ArgumentError hphc_fd(range(1, 512; step = 0.25), params, g)
            @test_throws ArgumentError hphc_fd(range(0, 512; step = -0.25), params, g)
            @test_throws ArgumentError hphc_fd(range(0, 0; length = 1), params, g)
        end

        @testset "parameter validation" begin
            @test_throws ArgumentError hphc_fd(f, merge(params, (; not_a_parameter = 1.0)), g)
        end

        @testset "output" begin
            h = hphc_fd(f, params, g)
            @test h.plus isa AbstractVector{<:Complex}
            @test h.cross isa AbstractVector{<:Complex}
            @test length(h.plus) == length(f)
            @test length(h.cross) == length(f)

            f_min = minimum_frequency(g)
            f_max = something(maximum_frequency(g), Float64(last(f)))
            outside = [!(f_min <= fi <= f_max) for fi in f]
            inside = .!outside
            @test any(outside) && any(inside)
            for x in (h.plus, h.cross)
                @test all(iszero, x[outside])
                @test any(!iszero, x[inside])
                @test all(isfinite, x)
            end
        end
    end
end
