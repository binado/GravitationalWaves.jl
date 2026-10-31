using TestItems

# Evaluated afresh in each test item, including mutable fixtures and dummy types.
@testsnippet Common begin
    using Test
    using GWSeries
    using DimensionalData
    using AbstractFFTs
    using FFTW
    using LinearAlgebra
    using Random
    # Storage traits and static arrays come from DimensionalData's own extensions.
    using ArrayInterface
    using StaticArrays

    include("prototypes.jl")
end
