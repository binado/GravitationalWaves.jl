# Exercise array inputs and every numerical public API with dynamic geometry.
function numerical_api(x, f, hp, hc)
    d = ad_detector(x)
    other = detectors.L1
    sky = (; ra = x[1], dec = x[2], psi = x[3], geocent_time = 1126259462 + x[4])
    return (
        detector_tensor(d), d.vertex, d.xarm, d.yarm,
        greenwich_mean_sidereal_time(sky.geocent_time),
        antenna_pattern(sky, d),
        antenna_pattern(x[1], x[2], x[3], x[5], d),
        time_delay_from_geocenter(sky, d),
        time_delay_from_geocenter(x[1], x[2], x[5], d),
        time_delay(sky, d, other), light_travel_time(d, other),
        project_fd(
            f, (; plus = hp, cross = hc), sky, d; start_time = 1126259462.0,
            reference_time = 1126259462 + x[5]
        ),
    )
end

function compare_host(actual, expected)
    if expected isa Union{Tuple, NamedTuple}
        foreach(compare_host, values(actual), values(expected))
    elseif expected isa AbstractArray
        # StaticArrays return concrete device scalars; device arrays materialize as numbers.
        @test Reactant.to_number.(Array(actual)) ≈ expected rtol = 2.0e-10 atol = 2.0e-10
    else
        @test Reactant.to_number(actual) ≈ expected rtol = 2.0e-10 atol = 2.0e-10
    end
    return nothing
end

@testset "reused array executable" begin
    samples = (AD_INPUT, AD_INPUT .+ [0.1, 0.02, -0.03, 0.125, 0.25, 0.01, -0.02, 0.04, -0.01, 0.001, 0.002, 10.0])
    compiled_api(x, f, hp, hc) = Reactant.@allowscalar numerical_api(x, f, hp, hc)
    args = Reactant.to_rarray.((AD_INPUT, AD_F, AD_POL.plus, AD_POL.cross))
    executable = Reactant.@compile compiled_api(args...)
    for sample in samples
        f = sample === first(samples) ? AD_F : [0.0, 19.5, 35.0, 110.0]
        hp, hc = AD_POL.plus .* sample[1], AD_POL.cross .* sample[2]
        args = Reactant.to_rarray.((sample, f, hp, hc))
        compare_host(executable(args...), numerical_api(sample, f, hp, hc))
    end
end

# Cartesian constructor and its geometry derivatives also remain dynamic.
function cartesian_loss(x)
    u = (cos(x[1]), sin(x[1]), zero(x[1]))
    v = (zero(x[2]), cos(x[2]), sin(x[2]))
    d = Detector((x[3], x[4], x[5]), u, v, (3000, 4000); check = false)
    D = detector_tensor(d)
    return D[1, 2] + D[2, 3] + time_delay_from_geocenter(1.0, -0.3, 2.0, d)
end

@testset "Cartesian geometry gradients" begin
    objective(x) = Reactant.@allowscalar cartesian_loss(x)
    samples_cartesian = ([0.3, 1.2, 1.0e6, 2.0e6, 3.0e6], [0.5, 1.0, 2.0e6, -1.0e6, 4.0e6])
    tracked = Reactant.to_rarray(first(samples_cartesian))
    executable = Reactant.@compile Enzyme.gradient(Enzyme.Reverse, objective, tracked)
    for sample in samples_cartesian
        g = Array(only(executable(Enzyme.Reverse, objective, Reactant.to_rarray(sample))))
        @test g ≈ ForwardDiff.gradient(cartesian_loss, sample) rtol = 1.0e-10 atol = 1.0e-11
    end
end

@testset "reused GMST executable across leap-table intervals" begin
    objective(x) = Reactant.@allowscalar greenwich_mean_sidereal_time(x[1])
    epochs = (
        0.25, 46828799.25, 46828801.25, 1119744015.25,
        1119744017.25, 1167264016.25, 1167264018.25, 1400000000.25,
    )
    tracked = Reactant.to_rarray([first(epochs)])
    value = Reactant.@compile objective(tracked)
    gradient = Reactant.@compile Enzyme.gradient(Enzyme.Reverse, objective, tracked)
    for epoch in epochs
        x = [epoch]
        args = Reactant.to_rarray(x)
        @test Reactant.to_number(value(args)) ≈ objective(x) atol = 3.0e-11 rtol = 0
        @test Array(only(gradient(Enzyme.Reverse, objective, args))) ≈
            ForwardDiff.gradient(objective, x) rtol = 1.0e-12
    end
end
