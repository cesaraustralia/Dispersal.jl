"Extends AbstractNeighborhood for for dispersal kernel neighborhoods"
abstract type AbstractDispersalKernel <: AbstractNeighborhood end

@mix @columns struct Kernel{N}
    neighborhood::N | DispersalKernel() | true  | _ | "Neighborhood to disperse to or from"
end

"""
A neighborhood built from a dispersal kernel function and a cell radius.
Can be built directly by passing in the array, radius and overflow
arguments, but preferably use the keyword constructor to build the array from
a dispersal kernel function.
$(FIELDDOCTABLE)
"""
@description @limits @flattenable struct DispersalKernel{F,P,K,C,I} <: AbstractDispersalKernel
    f::F        | false | _           | "Kernel function"
    param::P    | true  | (0.0, 10.0) | "Parameter for dispersal kernel function"
    kernel::K   | false | _           | "Kernal matrix"
    temp::K     | false | _           | "Temp matrix for current neighborhood"
    cellsize::C | false | (0.0, 10.0) | "Simulation cell size"
    radius::I   | false | (1, 10)     | "Kernel radius"
    function DispersalKernel{F,P,K,C,I}(f::F, param::P, init_kernel::K, temp_kernel::K, 
                    cellsize::C, radius::I) where {F,P,K,C,I} 
        kernel = build_dispersal_kernel(f, param, init_kernel, cellsize, radius)
        temp_kernel = similar(kernel)
        k = typeof(kernel)
        new{F,P,k,C,I}(f, param, kernel, temp_kernel, cellsize, radius)
    end
end

"Constructor that generates the temp kernel from the kernel"
DispersalKernel(f::F, param::P, kernel::K, cellsize::C, radius::I) where {F,P,K,C,I} =
    DispersalKernel{F,P,K,C,I}(f, param, kernel, similar(kernel), cellsize, radius)

"""
    DispersalKernel(; f=exponential, param=1.0, init=[], cellsize=1.0, radius=3)
Constructor for DispersalKernel, accepting keyword arguments.
"""
DispersalKernel(; f=exponential, param=1.0, init=[], cellsize=1.0, radius=3) = 
    DispersalKernel(f, param, init, cellsize, radius)

build_dispersal_kernel(f, params, init, cellsize, r) = begin
    params = typeof(params) <: Tuple ? params : (params,)
    sze = 2r + one(r)
    kernel = zeros(Float64, sze, sze)
    # Paper: l. 97
    for y = -r:r, x = -r:r
        d = sqrt(y^2 + x^2) * cellsize
        kernel[y+r+one(y), x+r+one(x)] = f(d, params...)
    end
    # Normalise
    kernel ./= sum(kernel)
    typeof(init).name.wrapper(kernel)
end

Cellular.radius(hood::DispersalKernel) = hood.radius
Cellular.temp_neighborhood(hood::DispersalKernel) = hood.temp

# Paper: l. 96 TODO make this a type with the parameter self contained?
exponential(d, a) = exp(-d / a)
