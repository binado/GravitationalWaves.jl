using TestItems

# Evaluated afresh in each test item, including mutable fixtures and dummy types.
@testsnippet Common begin
    using Test
    using GWMatchedFiltering
    using LinearAlgebra
    using Random
    using ForwardDiff

    include("whitening_helpers.jl")
end
