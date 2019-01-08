"Kernel neighbourhoods for dispersal"
abstract type AbstractDispersalKernel <: AbstractNeighborhood end

@mix @columns struct Kernel{N}
    "Neighborhood to disperse to or from"
    neighborhood::N = DispersalKernel() | true  | _
end

"""
A neighborhood built from a dispersal kernel function and a cell radius.
Can be built directly by passing in the array, radius and overflow
arguments, but preferably use the keyword constructor to build the array from
a dispersal kernel function.
"""
@limits @flattenable struct DispersalKernel{F,P,K,C,I} <: AbstractDispersalKernel
    f::F        | false | _
    param::P    | true  | (0.0, 10.0)
    kernel::K   | false | _
    temp::K     | false | _
    cellsize::C | false | (0.0, 10.0)
    radius::I   | false | (1, 10)
    function DispersalKernel{F,P,K,C,I}(f::F, param::P, init_kernel::K, temp_kernel::K, 
                    cellsize::C, radius::I) where {F,P,K,C,I} 
        kernel = build_dispersal_kernel(f, param, init_kernel, cellsize, radius)
        temp_kernel = similar(kernel)
        k = typeof(kernel)
        new{F,P,k,C,I}(f, param, kernel, temp_kernel, cellsize, radius)
    end
end

DispersalKernel(f::F, param::P, kernel::K, cellsize::C, radius::I) where {F,P,K,C,I} =
    DispersalKernel{F,P,K,C,I}(f, param, kernel, similar(kernel), cellsize, radius)

"""
    DispersalKernel(; dir=:inwards, f=exponential, param=1.0, init=[], cellsize=1.0, radius=Int64(3), overflow=Skip())
Constructor for neighborhoods, using a dispersal kernel function and a cell radius.

### Keyword Arguments:
- `f::Function`: any function that accepts a Number argument and returns a propbability between 0.0 and 1.0
- `radius::Integer`: a positive integer
- `overflow = Skip()
"""
DispersalKernel(; f=exponential, param=1.0, init=[], cellsize=1.0, radius=3) = 
    DispersalKernel(f, param, init, cellsize, radius)


Cellular.radius(hood::DispersalKernel) = hood.radius
Cellular.temp_neighborhood(hood::DispersalKernel) = hood.temp


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

# Paper: l. 96
exponential(d, a) = exp(-d / a)
