using CuArrays,
      GPUArrays,
      CUDAnative

import CUDAnative: cudaconvert

pressure(model, source::CuDeviceArray, cc, randomstate, args...) = begin
    rnd = spec_rand(source, Float64, randomstate)
    CUDAnative.pow(rnd, model.prob_threshold) > (one(cc) - cc) / one(cc)
end

cudafields(x::T) where T = cudaconvert.(getfield.(x, fieldnames(T)))

cudaconvert(x::SuitabilityLayer) = SuitabilityLayer(cudaconvert(x.data))
cudaconvert(x::HumanLayer) = HumanLayer(cudaconvert(x.data))
cudaconvert(x::SuitabilitySequence) = SuitabilitySequence(cudaconvert(x.data), x.timespan)
cudaconvert(x::DispersalNeighborhood{T}) where T = begin
    fields = (x.f, x.param, cudaconvert(CuArray(x.kernel)), x.cellsize, x.radius, x.overflow)
    DispersalNeighborhood{T, typeof.(fields)...}(fields...)
end

build_dispersal_kernel(f, params, init::CuDeviceArray, cellsize, r) = begin
    kernel = build_dispersal_kernel(f, params, [], cellsize, r) 
    cudaconvert(CuArray(kernel))
end


spec_rand(source::CuDeviceArray, typ, randomstate, args...) = GPUArrays.gpu_rand(typ, CuArrays.CuKernelState(), randomstate)
