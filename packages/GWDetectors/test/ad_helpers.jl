# Shared by the default suite and every heavy AD backend. Local epoch offsets keep
# finite differences resolvable, while GMST still receives full GPS epochs.
const AD_INPUT = [1.2, -0.4, 0.7, 0.25, 0.15, 0.6, -1.1, 0.2, 1.3, 0.03, -0.02, 50.0]
const AD_F = [0.0, 17.25, 31.0, 93.75]
const AD_POL = (;
    plus = [1.0 + 0.2im, -0.3 + 0.7im, 0.5 - 0.2im, -0.1 + 0.4im],
    cross = [0.3 - 0.8im, 0.9 + 0.1im, -0.4 + 0.6im, 0.7 + 0.2im],
)

function ad_detector(x)
    return Detector(
        "custom"; latitude = x[6], longitude = x[7],
        xarm_azimuth = x[8], yarm_azimuth = x[9], xarm_tilt = x[10],
        yarm_tilt = x[11], elevation = x[12], arm_lengths = (3000, 4000), check = false
    )
end

function geometry_loss(x)
    d = ad_detector(x)
    D = detector_tensor(d)
    return sum(d.vertex) / 6378137 + sum(d.xarm) + sum(d.yarm) + D[1, 2] + 2D[2, 3]
end

function antenna_loss(x)
    a = antenna_pattern(x[1], x[2], x[3], x[5], ad_detector(x))
    return a.plus + 0.3a.cross
end

function timing_loss(x)
    d = ad_detector(x)
    sky = (; ra = x[1], dec = x[2], geocent_time = 1126259462 + x[4])
    other = Detector("L1")
    return 10time_delay_from_geocenter(sky, d) +
        2time_delay(sky, d, other) + light_travel_time(d, other) +
        0.001greenwich_mean_sidereal_time(sky.geocent_time)
end

function projection_loss(x)
    sky = (; ra = x[1], dec = x[2], psi = x[3], geocent_time = 1126259462 + x[4])
    h = project_fd(
        AD_F, AD_POL, sky, ad_detector(x); start_time = 1126259462.0,
        reference_time = 1126259462 + x[5]
    )
    return sum(real.(h) .+ 0.3 .* imag.(h))
end
