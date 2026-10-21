if !isdefined(@__MODULE__, :time_prototype)
    include(joinpath(@__DIR__, "../../packages/GWSeries/test/prototypes.jl"))
end
const VECTOR_DOMAIN = frequency_prototype(; frequencies = 0:0.5:4, epoch = 1126259460.0, channels = (:H1,))
const VECTOR_METRIC = FrequencyDomainMetric(VECTOR_DOMAIN, ones(9); nsamples = 16)
const VECTOR_F = collect(frequencies(VECTOR_DOMAIN))
const VECTOR_DATA = complex.(sin.(1:9), cos.(1:9))

function vector_loss(x)
    signal = x[1] .* VECTOR_DATA .* exp.(1.0im .* x[2] .* VECTOR_F)
    return real(matched_filter_snr(signal, VECTOR_DATA, VECTOR_METRIC)) +
        log_likelihood_ratio(signal, VECTOR_DATA, VECTOR_METRIC) +
        residual_squared(VECTOR_DATA, signal, VECTOR_METRIC)
end

@testset "Single-channel vector kernels and fixed domains" begin
    check_data(VECTOR_DOMAIN, VECTOR_DATA)
    check_template(VECTOR_DATA, VECTOR_METRIC)
    sample = [0.8, 0.2]
    expected = ForwardDiff.gradient(vector_loss, sample)
    for backend in (AutoEnzyme(; function_annotation = Enzyme.Const), AutoMooncake(; config = nothing))
        @test DifferentiationInterface.gradient(vector_loss, backend, sample) ≈ expected rtol = 1.0e-7 atol = 1.0e-8
    end
    objective(x) = Reactant.@allowscalar vector_loss(x)
    tracked = Reactant.to_rarray(sample)
    executable = Reactant.@compile objective(tracked)
    gradient = Reactant.@compile Enzyme.gradient(Enzyme.Reverse, objective, tracked)
    @test Reactant.to_number(executable(tracked)) ≈ vector_loss(sample) rtol = 1.0e-8
    @test Array(only(gradient(Enzyme.Reverse, objective, tracked))) ≈ expected rtol = 1.0e-7 atol = 1.0e-8
end
