"Kernel neighbourhoods for dispersal"
abstract type AbstractDispersalKernel <: AbstractNeighborhood end

@mix struct Kernel{N}
    "Neighborhood to disperse to or from"
    neighborhood::N = DispersalKernel() | true  | _
end

"""
A neighborhood built from a dispersal kernel function and a cell radius.
Can be built directly by passing in the array, radius and overflow
arguments, but preferably use the keyword constructor to build the array from
a dispersal kernel function.
"""
@limits @flattenable struct DispersalKernel{F,P,K,C,I,O} <: AbstractDispersalKernel
    f::F        | false | _
    param::P    | true  | (0.0, 10.0)
    kernel::K   | false | _
    temp::K     | false | _
    cellsize::C | false | (0.0, 10.0)
    radius::I   | false | (1, 10)
    overflow::O | false | _
    function DispersalKernel{F,P,K,C,I,O}(f::F, param::P, init_kernel::K, 
                                                  cellsize::C, radius::I, overflow::O
                                                 ) where {T,F,P,K,C,I,O}
        kernel = build_dispersal_kernel(f, param, init_kernel, cellsize, radius)
        temp = similar(kernel)
        k = typeof(kernel)
        new{F,P,k,C,I,O}(f, param, kernel, temp, cellsize, radius, overflow)
    end
end

"""
    DispersalKernel(; dir=:inwards, f=exponential, param=1.0, init=[], cellsize=1.0, radius=Int64(3), overflow=Skip())
Constructor for neighborhoods, using a dispersal kernel function and a cell radius.

### Keyword Arguments:
- `f::Function`: any function that accepts a Number argument and returns a propbability between 0.0 and 1.0
- `radius::Integer`: a positive integer
- `overflow = Skip()
"""
DispersalKernel(; f=exponential, param=1.0, init=[], cellsize=1.0, 
                      radius=Int64(3), overflow=Skip()) = begin
    DispersalKernel{typeof.((f, param, init, cellsize, radius, overflow))...
                         }(f, param, init, cellsize, radius, overflow)
end


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
