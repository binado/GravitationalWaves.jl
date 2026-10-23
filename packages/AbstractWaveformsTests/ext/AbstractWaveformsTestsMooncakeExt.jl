module AbstractWaveformsTestsMooncakeExt

using AbstractWaveformsTests: FixedObjective
using GWMatchedFiltering: FrequencyDomainMetric
using Mooncake

# Captures of a FixedObjective are fixed configuration, never aliases of the differentiated input.
Mooncake.tangent_type(::Type{<:FixedObjective}) = Mooncake.NoTangent
# A prepared metric is fixed configuration, and Mooncake cannot copy the TwicePrecision
# frequency metadata stored inside it.
Mooncake.tangent_type(::Type{<:FrequencyDomainMetric}) = Mooncake.NoTangent

end
