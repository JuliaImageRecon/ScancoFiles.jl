module ScancoFiles

using ImageMetadata
using Unitful: cm

# see also https://github.com/KitwareMedical/ITKIOScanco/blob/master/src/itkScancoImageIO.cxx

export loadISQ, loadRSQ

function readScancoInt(fd)
  return read(fd, Int32)
end

function readScancoFloat32(fd)
  return read(fd, Float32)
end

function readScancoFloat64(fd)
  return read(fd, Float64)
end

function readScancoString(fd, len)
  return String( read!(fd, Array{UInt8}(undef, len)) )
end


function load(filename::AbstractString; kargs...)
  filenamebase, ext = splitext(filename)

  if lowercase(ext) == ".isq"
    return loadISQ(filename; kargs...)
  elseif lowercase(ext) == ".rsq"
    return loadRSQ(filename; kargs...)
  else
    error("cannot read $(ext) files")
  end
end

function loadISQ(filename::AbstractString)
  filenamebase, ext = splitext(filename)

  I = open(filename, "r")  do fd
    params = Dict{String, Any}()
    params["version"] = readScancoString(fd, 16)
    params["dataType"] = readScancoInt(fd)
    params["numBytes"] = readScancoInt(fd)
    params["numBlocks"] = readScancoInt(fd)
    params["patientIndex"] = readScancoInt(fd)
    params["scannerID"] = readScancoInt(fd)
    # TODO read date correctly
    read!(fd, Array{UInt8}(undef, 8))
    params["pixdim"] = ntuple(_->readScancoInt(fd), 3)
    params["physdim"] = ntuple(_->readScancoInt(fd)*1e-4cm, 3)

    isRAD = (params["dataType"] == 9 || params["physdim"][3] == 0)
    if isRAD
      error("cannot read rad files")
    end

    params["sliceThickness"] = readScancoInt(fd)*1e-3
    params["sliceIncrement"] = readScancoInt(fd)*1e-3
    params["startPosition"] = readScancoInt(fd)*1e-3
    # this->m_EndPosition = this->m_StartPosition + physdim[2] * 1e-3 * (pixdim[2] - 1) / pixdim[2];
    params["dataRange"] = ntuple(_->readScancoInt(fd), 2)
    params["muScaling"] = readScancoInt(fd)
    params["numberOfSamples"] = readScancoInt(fd)
    params["numberOfProjections"] = readScancoInt(fd)
    params["scanDistance"] = readScancoInt(fd)*1e-3
    params["scannerType"] = readScancoInt(fd)
    params["sampleTime"] = readScancoInt(fd)*1e-3
    params["measurementIndex"] = readScancoInt(fd)
    params["site"] = readScancoInt(fd)
    params["referenceLine"] = readScancoInt(fd)*1e-3
    params["reconstructionAlg"] = readScancoInt(fd)
    params["patientName"] = strip(readScancoString(fd, 40))
    params["energy"] = readScancoInt(fd)*1e-3
    params["intensity"] = readScancoInt(fd)*1e-3

    seek(fd, 508)

    params["dataOffset"] = readScancoInt(fd)

    params["headerSize"] = (params["dataOffset"] + 1) * 512

    

    dtype = Int16 # TODO: match this
  
    seek(fd,512)
    specialMode = read(fd,Int16) == 2573    

    if !specialMode
      seek(fd, params["headerSize"])
      I = read!(fd, Array{dtype}(undef, params["pixdim"]))
    else
      seek(fd,0)
      Q = read!(fd, Array{Int16}(undef, filesize(fd)÷2))
      Y = deleteat!(Q, 257:257:length(Q))
      j = params["headerSize"] ÷ 2 + 1
      I = collect(reshape( Y[j:j+prod(params["pixdim"])-1], params["pixdim"]))
    end

    # TODO: Scale I
    return ImageMeta(I, ImageMetadata.to_dict(params))
  end
  
  return I
end


function loadRSQ(filename::AbstractString; sinogramPreprocessing=true)
  filenamebase, ext = splitext(filename)
  
  I = open(filename, "r")  do fd
    params = Dict{String, Any}()
    params["version"] = readScancoString(fd, 16)
    params["dataType"] = readScancoInt(fd)
    params["numBytes"] = readScancoInt(fd)
    params["numBlocks"] = readScancoInt(fd)
    params["patientIndex"] = readScancoInt(fd)
    params["scannerID"] = readScancoInt(fd)
    # TODO read date correctly
    read!(fd, Array{UInt8}(undef, 8))
    params["pixdim"] = ntuple(_->readScancoInt(fd), 3)
    params["physdim"] = ntuple(_->readScancoInt(fd)*1e-4cm, 3)

    isRAD = (params["dataType"] ==9 || params["physdim"][3] == 0)
    if isRAD
      error("cannot read rad files")
    end

    params["sliceThickness"] = readScancoInt(fd)*1e-3
    params["sliceIncrement"] = readScancoInt(fd)*1e-3
    params["startPosition"] = readScancoInt(fd)*1e-3
    # this->m_EndPosition = this->m_StartPosition + physdim[2] * 1e-3 * (pixdim[2] - 1) / pixdim[2];
    params["dataRange"] = ntuple(_->readScancoInt(fd), 2)
    params["muScaling"] = readScancoInt(fd)
    params["numberOfSamples"] = readScancoInt(fd)
    params["numberOfProjections"] = readScancoInt(fd)
    params["detectorWidth"] = readScancoInt(fd)*1e-4cm
    params["detectorHeight"] = readScancoInt(fd)*1e-4cm
    params["centerPixel"] = ntuple(_->readScancoInt(fd)*1e-4cm, 2)
    params["detectorSourceDist"] = readScancoInt(fd)*1e-4cm
    params["detectorCenterDist"] = params["detectorSourceDist"] - readScancoInt(fd)*1e-4cm
    readScancoInt(fd)
    read!(fd, Array{UInt8}(undef, 12 * 4))
    params["dataPixelR"] = readScancoInt(fd)
    readScancoInt(fd)
    readScancoInt(fd)
    params["ioRecord"] = readScancoInt(fd)
    params["darkRecord"] = readScancoInt(fd)
    params["dataRecord"] = readScancoInt(fd)
    readScancoInt(fd)
    params["integrationTime"] = readScancoInt(fd)
    params["patientName"] = strip(readScancoString(fd, 40))
    params["energy"] = readScancoInt(fd)
    params["intensity"] = readScancoInt(fd)

    seek(fd, 480)
    params["orbitStartOffset"] = readScancoInt(fd)*1e-3

    seek(fd, 508)
    params["dataOffset"] = readScancoInt(fd)

    params["headerSize"] = (params["dataOffset"] + 1) * 512

    seek(fd, params["headerSize"])

    dtype = Int16 
  
    seek(fd,512)
    specialMode = read(fd,Int16) == 2573    

    if !specialMode
      seek(fd, params["headerSize"])
      I = read!(fd, Array{dtype}(undef, params["pixdim"]))
    else
      seek(fd,0)
      Q = read!(fd, Array{Int16}(undef, filesize(fd)÷2))
      Y = deleteat!(Q, 257:257:length(Q))
      j = params["headerSize"] ÷ 2 + 1
      I = collect(reshape( Y[j:j+prod(params["pixdim"])-1], params["pixdim"]))
    end
    
    # derived parameters
    params["detectorOffset"] = 0      
    d = params["detectorWidth"] / size(I,1) 
    a = sqrt((2*params["detectorSourceDist"])^2 - params["detectorWidth"]^2)*0.5
    β = 2*asin(params["detectorWidth"]/(2*a))/pi*180 
    params["fanAngle"] = β
    params["angleRange"] = 180+params["fanAngle"] #γ

    params["calibrationScans"] = I[:,1:2,:]
    I = I[:,3:end,:]
    reverse!(I, dims=2)

    if sinogramPreprocessing
      C = params["calibrationScans"]
      I_ = similar(I, Float32)
      for l=1:size(I,3)
        I_[:,:,l] = (-log.((I[:,:,l] .- C[:,1,l]) ./ (C[:,2,l] .- C[:,1,l]))) * 2^13 # magic factor
      end
      I = I_
    end

    return ImageMeta(I, ImageMetadata.to_dict(params))
  end
 
  return I
end


end
