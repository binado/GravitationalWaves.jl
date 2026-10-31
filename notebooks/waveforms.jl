### A Pluto.jl notebook ###
# v1.1.0

using Markdown
using InteractiveUtils

# ╔═╡ 399eec00-d545-11f1-bc3b-1ba75567c17e
begin
	import Pkg
	Pkg.activate(@__DIR__)

	using CairoMakie
	using FrequencyDomainWaveforms
	using GWSeries
end

# ╔═╡ 28a57670-87d2-417a-8336-27ea28b1ec7e
md"""
# Generating waveforms
"""

# ╔═╡ f372e778-131c-45e5-bc30-9f5d450f6653
params = (;
	mass_1 = 30.0, # detector-frame solar masses
	mass_2 = 30.0,
	luminosity_distance = 1000.0, # Mpc
	theta_jn = 0.4, # inclination in radians
	phase = 0.0, # coalescence phase in radians
)

# ╔═╡ 5ccaa5c5-5ef2-40ff-87fc-486de091927f
waveform_generator = IMRPhenomXAS(; minimum_frequency = 20.0, maximum_frequency = 512.0)

# ╔═╡ 9572a2c4-79c4-4aa3-920e-2baa7595fb48
frequency_grid = 0.0:0.25:512.0

# ╔═╡ 493f4649-c69d-4883-9d4c-2e46272ad545
prototype = FrequencySeries(zeros(ComplexF64, length(frequency_grid));
	epoch = 0.0, frequency_spacing = step(frequency_grid), channels = (:prototype,))

# ╔═╡ 61d813ee-20a9-4c5a-ae2e-85a23bd0cf01
begin
	validate_parameters(frequency_grid, params, waveform_generator)
	waveform = hphc_fd(prototype, params, waveform_generator)
end

# ╔═╡ 46a7670e-2380-4b30-99f6-041dfe14e123
lines(waveform; representation = :magnitude, axis = (; title = "IMRPhenomXAS polarizations", ylabel = "Strain magnitude (1/Hz)"))

# ╔═╡ Cell order:
# ╠═28a57670-87d2-417a-8336-27ea28b1ec7e
# ╠═399eec00-d545-11f1-bc3b-1ba75567c17e
# ╠═f372e778-131c-45e5-bc30-9f5d450f6653
# ╠═5ccaa5c5-5ef2-40ff-87fc-486de091927f
# ╠═9572a2c4-79c4-4aa3-920e-2baa7595fb48
# ╠═493f4649-c69d-4883-9d4c-2e46272ad545
# ╠═61d813ee-20a9-4c5a-ae2e-85a23bd0cf01
# ╠═46a7670e-2380-4b30-99f6-041dfe14e123
