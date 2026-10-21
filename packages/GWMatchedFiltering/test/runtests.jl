using Test
using GWSeries
using GWMatchedFiltering
using LinearAlgebra
using Random
using ForwardDiff

include("test_metric.jl")
include("test_statistics.jl")
include("ad_helpers.jl")
include("test_ad.jl")
include("test_reference.jl")
include("test_domains.jl")
include("test_quality.jl")
