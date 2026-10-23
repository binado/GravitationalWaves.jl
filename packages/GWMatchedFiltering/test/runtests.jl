using Test
using GWSeries
using GWMatchedFiltering
using LinearAlgebra
using Random
using ForwardDiff

include("whitening_helpers.jl")
include("test_metric.jl")
include("test_statistics.jl")
include("test_whitened.jl")
include("ad_helpers.jl")
include("test_ad.jl")
include("test_reference.jl")
include("test_series.jl")
include("test_quality.jl")
