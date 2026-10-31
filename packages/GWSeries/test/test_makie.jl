@testitem "Makie series plots" begin
    using GWSeries, Makie, Test

    x = GWSeries.TimeSeries([1.0 2.0; 3.0 4.0]; epoch = 1126259460.125, sampling_interval = 0.25, channels = (:H1, :L1))
    fig, ax, p = lines(x)
    @test length(ax.scene.plots) == 2
    @test ax.xlabel[] == "Time since GPS 1.126259460125e9 (s)"
    @test ax.scene.plots[1].label[] == "H1"
    @test ax.scene.plots[2].label[] == "L1"
    @test p[1][][2][1] == 0.25
    ax2 = Axis(Figure()[1, 1])
    lines!(ax2, x)
    @test length(ax2.scene.plots) == 2

    h = FrequencySeries(reshape(ComplexF64[3 + 4im, -1 + 2im], 2, 1); epoch = 0, frequency_spacing = 0.5, channels = (:H1,))
    @test_throws ArgumentError lines(h)
    @test_throws ArgumentError lines(h; representation = :unknown)
    for (representation, expected) in ((:magnitude, 5.0), (:phase, angle(3 + 4im)), (:real, 3.0), (:imaginary, 4.0))
        _, axis, plot = lines(h; representation)
        @test axis.xlabel[] == "Frequency (Hz)"
        @test plot[1][][1][2] ≈ expected
        @test plot[1][][2][1] == 0.5
    end
    real_h = FrequencySeries([-2.0, 3.0]; epoch = 0, frequency_spacing = 0.5, channels = (:H1,))
    _, _, real_plot = lines(real_h)
    @test real_plot[1][][1][2] == -2
end
