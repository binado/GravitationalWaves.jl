using Test

# Peak-normalized error avoids division by tiny individual samples. Zero outputs
# are checked exactly instead of normalizing by zero.
function check_oracle(actual, expected; rtol)
    @test size(actual) == size(expected)
    @test all(isfinite, actual)
    @test all(isfinite, expected)
    peak = maximum(abs, expected)
    @test iszero(peak) ? all(iszero, actual) : maximum(abs, actual .- expected) <= rtol * peak
    return nothing
end
