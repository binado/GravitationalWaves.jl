using TestItems

@testitem "kernels receive plain ranges from series prototypes" setup = [Common] begin
    probe = GridProbe(Ref{Any}(nothing))
    fd = FrequencySeries(zeros(9, 1); epoch = 1126259460.0, frequency_spacing = 0.5, channels = (:H1,))
    hphc_fd(fd, (;), probe)
    @test probe.seen[] isa StepRangeLen
    @test probe.seen[] == frequencies(fd)
    td = TimeSeries(zeros(8, 1); epoch = 1126259460.125, sampling_interval = 0.25, channels = (:H1,))
    hphc_td(td, (;), probe)
    @test probe.seen[] isa StepRangeLen
    @test probe.seen[] == times(td)
end

@testitem "Aqua" setup = [Common, Setup_modes] begin
    Aqua.test_all(AbstractWaveforms)
end

@testitem "JET" setup = [Common, Setup_modes] begin
    JET.test_package(AbstractWaveforms; target_modules = (AbstractWaveforms,))
end

@testitem "capability fallbacks" setup = [Common] begin
    g = Dummy(nothing)
    @test spin_support(g) === :aligned
    @test supports_mode_array(g) === false
    @test is_differentiable(g) === false
end

@testitem "check_frequency_grid" setup = [Common] begin
    @test check_frequency_grid(0:0.5:10) == 0:0.5:10
    @test_throws ArgumentError check_frequency_grid(1:0.5:10)
    @test_throws ArgumentError check_frequency_grid(0:-0.5:-10)
    @test_throws ArgumentError check_frequency_grid(range(0, 0; length = 1))
end

@testitem "frequency_bounds" setup = [Common] begin
    f = 0:0.5:100
    @test frequency_bounds(Dummy(nothing), f) == (20.0, 100.0)
    @test frequency_bounds(Dummy(60.0), f) == (20.0, 60.0)
end

@testitem "show" setup = [Common] begin
    @test sprint(show, Dummy(60.0)) ==
        "Dummy(\"Dummy\"; reference_frequency=50.0, minimum_frequency=20.0, maximum_frequency=60.0)"
end

@testitem "quadrupole generator contract" setup = [Common] begin
    test_generator_contract(DummyQuadrupole(); params = (; theta_jn = 0.4), f = 0:0.5:100)
end

@testitem "unbounded generator contract" setup = [Common] begin
    test_generator_contract(Dummy(nothing); params = (mass_1 = 30.0, mass_2 = 30.0), f = 0:0.5:100)
end

@testitem "bounded generator contract" setup = [Common] begin
    test_generator_contract(Dummy(60.0); params = (mass_1 = 30.0, mass_2 = 30.0), f = 0:0.5:100)
end
