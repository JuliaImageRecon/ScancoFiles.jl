using ScancoFiles
using Test
using Scratch
using Unitful: cm

const datadir  = @get_scratch!("ScancoFiles")
@info "If you want to check the output of the tests, please head to $datadir."

filenameISQ = joinpath(datadir, "C0004255.ISQ")
if !isfile(filenameISQ)
  download("https://data.kitware.com/api/v1/file/591e56178d777f16d01e0d20/download", filenameISQ)
end
@info filenameISQ

@testset "ScancoFiles.jl" begin
    img = loadISQ(filenameISQ)

    @test img.version == "CTDATA-HEADER_V1"
    @test size(img) ==  (1024, 1024, 51)
    @test img.scannerID == 2135
    @test img.numberOfProjections == 500
    @test img.scannerType == 10
    @test img.intensity == 0.177
    @test img.dataRange == (-2813, 32767)
    @test img.pixdim == (1024, 1024, 51)
    @test img.reconstructionAlg == 3
    @test img.sliceThickness == 0.036000000000000004
    @test img.sliceIncrement == 0.036000000000000004
    @test img.measurementIndex == 4937
    @test img.referenceLine == 0.0
    @test img.site == 5
    @test img.dataOffset == 6
    @test img.patientName == "COLE-BPBP"                               
    @test img.numberOfSamples == 1024
    @test img.sampleTime == 400.0
    @test img.muScaling == 4096
    @test all(img.physdim .≈ (3.6864, 3.6864, 0.1836).*cm)
    @test img.numBlocks == 208903
    @test img.energy == 45.0
    @test img.dataType == 3
    @test img.headerSize == 3584
    @test img.scanDistance == 36.864000000000004
    @test img.startPosition == 75.0
    @test img.numBytes == 106958336
    @test img.patientIndex == 78

end
