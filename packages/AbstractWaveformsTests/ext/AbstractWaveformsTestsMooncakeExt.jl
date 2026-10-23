module AbstractWaveformsTestsMooncakeExt

using AbstractWaveformsTests: FixedObjective
using Mooncake

# Captures of a FixedObjective are fixed configuration, never aliases of the differentiated input.
Mooncake.tangent_type(::Type{<:FixedObjective}) = Mooncake.NoTangent

end
