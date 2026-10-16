import BackgroundCosmology
using BackgroundCosmology: cosmology, comoving_volume_element
using Distributions
using ForwardDiff
using Random

const _RS_PARAMS = (H0 = 67.66, Ωm = 0.30966, γ = 2.7, κ = 3.0, z_peak = 2.0, local_merger_rate = 17.0)

@testset "madau_dickinson_rate" begin
    γ, κ, z_peak = 2.7, 3.0, 2.0
    @test madau_dickinson_rate(0.0, γ, κ, z_peak) ≈ 1.0
    @test madau_dickinson_rate(0.0, γ, κ, z_peak, 17.0) ≈ 17.0
    for z in (0.5, 2.0, 3.0)
        by_hand = (1 + (1 + z_peak)^(-(γ + κ))) * (1 + z)^γ /
            (1 + ((1 + z) / (1 + z_peak))^(γ + κ))
        @test madau_dickinson_rate(z, γ, κ, z_peak) ≈ by_hand
    end
end

@testset "RedshiftDistribution table" begin
    d = madau_dickinson_redshift_distribution(; params = _RS_PARAMS, maximum_redshift = 2.0, n_grid = 101)
    c = cosmology(_RS_PARAMS)
    z = collect(range(0.0, 2.0; length = 101))
    @test redshift_grid(d) == z
    @test n_grid(d) == 101
    @test minimum(d) == 0.0 && maximum(d) == 2.0

    # Hand-built table: ψ(z)/(1+z) · dV_c/dz with dV_c/dz = 4π · (c/H₀) d_c² / E.
    dvc = @. 4π * comoving_volume_element(c, z)
    ψ = madau_dickinson_rate.(z, 2.7, 3.0, 2.0, 17.0)
    y = @. ψ / (1 + z) * dvc
    @test d.dist.y ≈ y rtol = 1.0e-9
    @test normalizer(d) === GWDistributions._trapz(d.dist.y, z)
    @test total_merger_rate(d) ≈ 1.0e-9 * GWDistributions._trapz(y, z) / JULIAN_YEAR_SEC rtol = 1.0e-9

    # Quadrature version over the same nodes agrees up to the trapezoid weights.
    nodes = collect(range(0.0, 2.0; length = 4001))
    weights = fill(nodes[2] - nodes[1], length(nodes))
    weights[[1, end]] ./= 2
    quad = total_merger_rate(nodes, weights, madau_dickinson_rate.(nodes, 2.7, 3.0, 2.0, 17.0), c)
    @test quad ≈ total_merger_rate(d) rtol = 1.0e-4

    @test merger_rate(d, 1.3) == madau_dickinson_rate(1.3, 2.7, 3.0, 2.0, 17.0)
    @test source_frame_distribution(d, 1.3, _RS_PARAMS) == merger_rate(d, 1.3)
    @test MadauDickinsonRedshiftDistribution === madau_dickinson_redshift_distribution

    # Support behaviour.
    @test isfinite(logpdf(d, 0.5))
    @test logpdf(d, -0.1) == -Inf
    @test logpdf(d, 2.1) == -Inf
    @test cdf(d, quantile(d, 0.3)) ≈ 0.3 atol = 1.0e-12
    samples = rand(MersenneTwister(1234), d, 100)
    @test all(s -> insupport(d, s), samples)
end

@testset "RedshiftDistribution distance accessors" begin
    d = madau_dickinson_redshift_distribution(; params = _RS_PARAMS, maximum_redshift = 2.0, n_grid = 101)
    c = cosmology(_RS_PARAMS)
    # On a node: exact; between nodes: linear; outside: clamped.
    z_node = redshift_grid(d)[41]
    @test luminosity_distance(d, z_node) ≈ BackgroundCosmology.luminosity_distance(c, z_node) rtol = 1.0e-8
    @test differential_comoving_volume(d, z_node) ≈ 4π * comoving_volume_element(c, z_node) rtol = 1.0e-8
    z_mid = 0.5 * (redshift_grid(d)[41] + redshift_grid(d)[42])
    @test luminosity_distance(d, z_mid) ≈
        (d.luminosity_distance_grid[41] + d.luminosity_distance_grid[42]) / 2
    @test luminosity_distance(d, 5.0) == last(d.luminosity_distance_grid)
    @test luminosity_distance(d, -1.0) == 0.0

    # Without xi_0 the GW distance is the EM distance.
    @test distance_ratio(d, 1.0) == 1.0
    @test gw_luminosity_distance(d, 1.0) == luminosity_distance(d, 1.0)

    p = merge(_RS_PARAMS, (xi_0 = 0.8, xi_n = 1.5))
    dm = madau_dickinson_redshift_distribution(; params = p, maximum_redshift = 2.0, n_grid = 101)
    expected = 0.8 + (1 - 0.8) * (1 + 1.0)^(-1.5)
    @test distance_ratio(dm, 1.0) ≈ expected
    @test gw_luminosity_distance(dm, 1.0) ≈ luminosity_distance(dm, 1.0) * expected
    # xi_0 = 1 is GR propagation.
    dgr = madau_dickinson_redshift_distribution(;
        params = merge(_RS_PARAMS, (xi_0 = 1.0, xi_n = 2.0)), maximum_redshift = 2.0, n_grid = 11
    )
    @test distance_ratio(dgr, 1.7) ≈ 1.0
end

@testset "RedshiftDistribution with w0waCDM" begin
    p = merge(_RS_PARAMS, (w0 = -0.9, wa = 0.1))
    d = madau_dickinson_redshift_distribution(; params = p, maximum_redshift = 2.0, n_grid = 51)
    @test isfinite(normalizer(d))
    @test isfinite(logpdf(d, 1.0))
end

@testset "RedshiftDistribution preserves AD" begin
    f = Ωm -> begin
        p = (H0 = 67.66, Ωm = Ωm, γ = 2.7, κ = 3.0, z_peak = 2.0, local_merger_rate = 17.0)
        d = madau_dickinson_redshift_distribution(; params = p, maximum_redshift = 2.0, n_grid = 51)
        total_merger_rate(d)
    end
    g = ForwardDiff.derivative(f, 0.31)
    @test isfinite(g) && g != 0.0

    h = γ -> begin
        p = merge(_RS_PARAMS, (γ = γ,))
        d = madau_dickinson_redshift_distribution(; params = p, maximum_redshift = 2.0, n_grid = 51)
        logpdf(d, 0.7)
    end
    @test isfinite(ForwardDiff.derivative(h, 2.7))
    grad = ForwardDiff.gradient(
        v -> begin
            p = merge(_RS_PARAMS, (Ωm = v[1], γ = v[2]))
            d = madau_dickinson_redshift_distribution(; params = p, maximum_redshift = 2.0, n_grid = 51)
            logpdf(d, 0.7) + quantile(d, 0.4)
        end, [0.31, 2.7]
    )
    @test all(isfinite, grad)
end
