using Distributions

# Values printed by astrogwb (JAX, x64) with H0 = 67.66, Ω_m = 0.30966, γ = 2.7, κ = 3.0,
# z_peak = 2.0, local_merger_rate = 17 on the default grid (z ∈ [0, 10], 1000 nodes).
@testset "astrogwb reference values" begin
    params = (H0 = 67.66, Ωm = 0.30966, γ = 2.7, κ = 3.0, z_peak = 2.0, local_merger_rate = 17.0)
    zs = [0.1, 0.5, 1.0, 2.5, 7.0]
    qs = [0.1, 0.5, 0.9]

    md = madau_dickinson_redshift_distribution(; params)
    @test logpdf.(md, zs) ≈ [
        -5.929663624193438, -2.6362611656691297, -1.409797634212804,
        -1.203173099105441, -4.6025482109387506,
    ] rtol = 1.0e-10
    @test total_merger_rate(md) ≈ 0.0023610007237638503 rtol = 1.0e-10
    @test normalizer(md) ≈ 74507516440250.08 rtol = 1.0e-10
    @test quantile.(md, qs) ≈ [1.042991124201453, 2.1384477336917347, 4.194063103626501] rtol = 1.0e-10
    @test luminosity_distance.(md, zs) ≈ [
        475.87930086060527, 2921.1747037917157, 6797.445511453407,
        20915.189653757578, 70627.64881799121,
    ] rtol = 1.0e-10

    # BackgroundCosmology's lookback time is Gauss–Legendre rather than closed form.
    td = madau_dickinson_time_delayed_redshift_distribution(;
        params, time_delay_distribution = PowerLawDelay(-1.0, 0.02, 13.0)
    )
    @test logpdf.(td, zs) ≈ [
        -4.982381442238876, -1.97868886371515, -1.0248566108918993,
        -1.36212533650176, -5.282544909878776,
    ] rtol = 1.0e-8
    @test total_merger_rate(td) ≈ 0.000834697325064233 rtol = 1.0e-8
    @test merger_rate.(td, zs) ≈ [
        20.019460801302156, 34.07727904834067, 52.31356969970644,
        44.381458683505095, 3.1123849424402463,
    ] rtol = 1.0e-8
    @test quantile.(td, qs) ≈ [0.843353070165012, 1.8205159957675763, 3.566552908303526] rtol = 1.0e-8

    t = [0.05, 1.0, 10.0]
    for (α, cdfs, icdfs, logps) in (
            (
                -0.5,
                [0.023724699543317438, 0.2478482804577828, 0.8720389749716531],
                [0.23798235124467007, 3.5099509756796388, 10.621982351244663],
                [-0.4377425393302501, -1.9356086761072455, -3.0869012226042685],
            ),
            (
                -1.0,
                [0.1414689889869091, 0.6039894546782979, 0.9594927614681934],
                [0.03822269731997362, 0.5099019513592785, 6.802241030334998],
                [1.1274791003239197, -1.8682531732300711, -4.170838266224117],
            ),
        )
        d = PowerLawDelay(α, 0.02, 13.0)
        @test cdf.(d, t) ≈ cdfs rtol = 1.0e-10
        @test quantile.(d, qs) ≈ icdfs rtol = 1.0e-10
        @test logpdf.(d, t) ≈ logps rtol = 1.0e-10
    end

    n2 = MaxOfTwoNormals(1.4, 0.3)
    v = [1.2, 1.5, 1.8]
    @test cdf.(n2, v) ≈ [0.0637524815168843, 0.39760422347177043, 0.825897047152145] rtol = 1.0e-10
    @test quantile.(n2, qs) ≈ [1.2565179402871511, 1.563485640685208, 1.8896656368850597] rtol = 1.0e-10
    @test logpdf.(n2, v) ≈ [-0.620414355514084, 0.46147680520454004, -0.006350013948201533] rtol = 1.0e-10

    u2 = MaxOfTwoUniforms(1.1, 1.4)
    v = [1.2, 1.8, 2.4]
    @test cdf.(u2, v) ≈ [0.005102040816326518, 0.25, 0.8622448979591836] rtol = 1.0e-10
    @test quantile.(u2, qs) ≈ [1.5427188724235732, 2.089949493661167, 2.4281566172707194] rtol = 1.0e-10
    @test logpdf.(u2, v) ≈ [-2.2823823856765277, -0.33647223662121295, 0.28256697178501045] rtol = 1.0e-10

    uc = UniformCosine()
    θ = [0.3, 1.5, 2.8]
    @test cdf.(uc, θ) ≈ [0.02233175543719701, 0.46463139916614854, 0.971111170334329] rtol = 1.0e-10
    @test quantile.(uc, qs) ≈ [0.6435011087932844, 1.5707963267948966, 2.498091544796509] rtol = 1.0e-10
    @test logpdf.(uc, θ) ≈ [-1.912165243777173, -0.6956553367515417, -1.7868073010115837] rtol = 1.0e-10
end
