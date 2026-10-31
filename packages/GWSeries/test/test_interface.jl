using Interfaces
@testset "DimensionalData array interface" begin
    x = TimeSeries(reshape(collect(1.0:16.0), 8, 2); epoch = 100, sampling_interval = 0.25, channels = (:H1, :L1))
    h = rfft(x)
    @test Interfaces.test(DimensionalData.DimArrayInterface, TimeSeries, [x])
    @test Interfaces.test(DimensionalData.DimArrayInterface, FrequencySeries, [h])
end
