using CuArrays,
      GPUArrays,
      CUDAnative

import CUDAnative: cudaconvert
import Cellular: broadcast_rule!

pressure(model, source::CuDeviceArray, cc, randomstate, args...) = begin
    rnd = spec_rand(source, Float64, randomstate)
    CUDAnative.pow(rnd, model.prob_threshold) > (one(cc) - cc) / one(cc)
end

cudafields(x::T) where T = cudaconvert.(getfield.(x, fieldnames(T)))

cudaconvert(x::T) where T <: AbstractLayer = T.name.wrapper(cudafields(x)...)
cudaconvert(x::T) where T <: AbstractSequence = T.name.wrapper(cudafields(x)...)
cudaconvert(x::SuitabilitySequence) = SuitabilitySequence(cudaconvert(x.data), x.timespan)
cudaconvert(x::T) where T <: AbstractModel = T.name.wrapper(cudafields(x)...)
cudaconvert(x::DispersalNeighborhood{T, F, P, K, C, I, O}) where {T, F, P, K, C, I, O} = begin
    fields = (x.f, x.param, cudaconvert(CuArray(x.kernel)), x.cellsize, x.radius, x.overflow)
    DispersalNeighborhood{T, typeof.(fields)...}(fields...)
end

build_dispersal_kernel(f, params, init::CuDeviceArray, cellsize, r) = begin
    kernel = build_dispersal_kernel(f, params, [], cellsize, r) 
    cudaconvert(CuArray(kernel))
end


spec_rand(source::CuDeviceArray, typ, randomstate, args...) = GPUArrays.gpu_rand(typ, CuArrays.CuKernelState(), randomstate)


build_cell_pop_index(i, j, ii, jj, human::CuDeviceArray, dist::CuDeviceArray) = begin 
    CUDAnative.exp(-dist[abs(i - ii) + 1, abs(j - jj) + 1]) * human[i, j]
end

broadcast_rule!(model::AbstractPartialModel, data::FrameData{A}, indices, args...
               ) where A <: CuDeviceArray = begin
    # Initialise the dest array
    data.dest .= data.source
    @cuda threads=100 cuda_rule!(Ref(model), Ref(data), data.source, indices, args...)
end

cuda_rule!(model, data, args...) = begin
    x = (blockIdx().x-1) * blockDim().x + threadIdx().x
    y = (blockIdx().y-1) * blockDim().y + threadIdx().y
    rule!(model, data, data.source[y, x], (y, x), args...)
end
