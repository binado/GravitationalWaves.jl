# Mirrors of the LAL C structs we touch. Layouts follow LALAtomicDatatypes.h and
# LALDatatypes.h; `test_fd.jl` asserts sizes and offsets.

struct LIGOTimeGPS
    gpsSeconds::Int32
    gpsNanoSeconds::Int32
end

# LALNumUnits == 7
struct LALUnit
    powerOfTen::Int16
    unitNumerator::NTuple{7, Int16}
    unitDenominatorMinusOne::NTuple{7, UInt16}
end

struct COMPLEX16Sequence
    length::UInt32
    data::Ptr{ComplexF64}
end

# LALNameLength == 64
struct COMPLEX16FrequencySeries
    name::NTuple{64, UInt8}
    epoch::LIGOTimeGPS
    f0::Float64
    deltaF::Float64
    sampleUnits::LALUnit
    data::Ptr{COMPLEX16Sequence}
end
