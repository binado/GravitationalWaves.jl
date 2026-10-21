using AbstractWaveforms: mode_sum, sylm_m2

struct SyntheticModes{T, F} <: HigherModesWaveformGenerator
    td::T
    fd::F
    calls::Base.RefValue{Int}
end
SyntheticModes(td, fd) = SyntheticModes(td, fd, Ref(0))
function AbstractWaveforms.waveform_modes_td(g::SyntheticModes, p; kwargs...)
    g.calls[] += 1
    return g.td
end
function AbstractWaveforms.waveform_modes_fd(f, p, g::SyntheticModes)
    g.calls[] += 1
    return g.fd
end
AbstractWaveforms.mode_angles(p, ::SyntheticModes) = (p.theta_jn, get(p, :azimuth, 0.0))
AbstractWaveforms.minimum_frequency(::SyntheticModes) = 0.0
AbstractWaveforms.maximum_frequency(::SyntheticModes) = nothing

struct QuadrupoleModes <: HigherModesWaveformGenerator end
AbstractWaveforms.minimum_frequency(::QuadrupoleModes) = 20.0
AbstractWaveforms.maximum_frequency(::QuadrupoleModes) = 60.0
function AbstractWaveforms.waveform_modes_fd(f, p, ::QuadrupoleModes)
    a = amplitude_phase_fd(f, p, DummyQuadrupole())
    h = a.amplitude .* cis.(a.phase)
    # Each branch contains one helicity; sqrt(16π/5) makes the usual face-on strain.
    z = zero.(h)
    scale = sqrt(16π / 5)
    return (; mode_indices = ((2, 2), (2, -2)), positive = hcat(z, scale .* h), negative = hcat(scale .* conj.(h), z))
end

# Independent O(N²) DFT; no FFT dependency or assumptions about signed-mode symmetry.
dft(x, bins) = [sum(x[n + 1] * cis(-2π * k * n / length(x)) for n in 0:(length(x) - 1)) for k in bins]

@testset "higher modes" begin
    @test supports_mode_array(QuadrupoleModes()) === true
    @test supports_mode_array(DummyQuadrupole()) === false

    @testset "harmonics and Python fixtures" begin
        for line in eachline(joinpath(@__DIR__, "harmonics.csv"))
            startswith(line, "#") && continue
            l, m, theta, phi, re, im = parse.(Float64, split(line, ','))
            # gwtools evaluates via cos(theta), losing near-pole precision.
            atol = min(theta, abs(π - theta)) < 1.0e-3 ? 2.0e-12 : 4.0e-15
            @test sylm_m2(Int(l), Int(m), theta, phi) ≈ complex(re, im) atol = atol rtol = 2.0e-13
        end
        for theta in (0.0, 1.0e-10, 0.8, π - 1.0e-10, Float64(π))
            for l in 2:5
                @test sum(abs2(sylm_m2(l, m, theta, 0.3)) for m in (-l):l) ≈ (2l + 1) / (4π)
            end
        end
        @test sylm_m2(2, 2, 0.2f0, 0.3f0) isa ComplexF32
        @test sylm_m2(2, 2, big"0.2", big"0.3") isa Complex{BigFloat}
    end

    labels = ((2, -1), (3, 0), (5, 3), (2, 2))
    n = 12
    modes = [complex(sin(0.7k + j), cos(0.2k - 2j)) for k in 1:n, j in eachindex(labels)]
    t = collect(0.0:(n - 1))
    td = (; t, mode_indices = labels, modes)
    p = (; theta_jn = 0.9, azimuth = 0.37)
    g = SyntheticModes(td, nothing)
    test_modes_td_contract(g; params = p)
    g.calls[] = 0
    out = hphc_td(g, p)
    @test g.calls[] == 1
    direct = [sum(sylm_m2(l, m, p.theta_jn, p.azimuth) * modes[k, j] for (j, (l, m)) in enumerate(labels)) for k in 1:n]
    @test out.t === t
    @test out.hplus ≈ real.(direct)
    @test out.hcross ≈ -imag.(direct)
    order = (3, 1, 4, 2)
    @test mode_sum(modes[:, collect(order)], labels[collect(order)] |> Tuple, 0.9, p.azimuth) ≈ direct

    f = 0.0:1.0:(n ÷ 2)
    positive = hcat((dft(modes[:, j], f) for j in eachindex(labels))...)
    negative = hcat((dft(modes[:, j], -f) for j in eachindex(labels))...)
    fd = (; mode_indices = labels, positive, negative)
    gfd = SyntheticModes(td, fd)
    outfd = hphc_fd(f, p, gfd)
    @test gfd.calls[] == 1
    reordered = (; mode_indices = Tuple(labels[collect(order)]), positive = positive[:, collect(order)], negative = negative[:, collect(order)])
    reordered_out = hphc_fd(f, p, SyntheticModes(td, reordered))
    @test outfd.plus ≈ reordered_out.plus
    @test outfd.cross ≈ reordered_out.cross
    @test isreal(outfd.plus[1]) && isreal(outfd.cross[1])
    @test positive[1, :] == negative[1, :]
    @test outfd.plus ≈ dft(out.hplus, f) atol = 1.0e-13
    @test outfd.cross ≈ dft(out.hcross, f) atol = 1.0e-13

    @testset "quadrupole reconstruction and masking" begin
        f = 0:0.5:100
        for theta_jn in (0.0, 0.4, π / 2, Float64(π)), phase in (0.0, 0.37, 1.2)
            p = (; theta_jn, phase)
            r = hphc_fd(f, p, QuadrupoleModes())
            a = amplitude_phase_fd(f, p, DummyQuadrupole())
            h = a.amplitude .* cis.(a.phase)
            quadrupole = hphc_fd(f, p, DummyQuadrupole())
            @test r.plus ≈ quadrupole.plus
            @test r.cross ≈ quadrupole.cross atol = 1.0e-14 rtol = 1.0e-13
            @test r.plus ≈ ((1 + cos(theta_jn)^2) / 2) .* h
            @test r.cross ≈ (-im * cos(theta_jn)) .* h atol = 1.0e-14 rtol = 1.0e-13
            test_modes_fd_contract(QuadrupoleModes(); params = p, f)
            quadrupole_modes = waveform_modes_fd(f, p, QuadrupoleModes())
            face = waveform_modes_fd(f, merge(p, (; theta_jn = 0.0)), QuadrupoleModes())
            @test quadrupole_modes == face
            for branch in (quadrupole_modes.positive, quadrupole_modes.negative)
                @test all(iszero, branch[(f .< 20) .| (f .> 60), :])
            end
        end
    end
end
