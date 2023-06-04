"""
    DispersalKernel <: AbstractKernelStencil

    DispersalKernel(; kw...)

Dispersal kernel for taking the dot product of the stencil and a matching
kernel of weights. May hold any `Stencil` object: the kernel
will be built to match the shape, using the `folumation`, `cellsize`
and `distancemethod`.

# Keyword Arguments

- `stencil`: Any DynamicGrids.jl `Stencil`, or an
    already constructed [`DispersalKernel`](@ref). Using this keyword means `radius` is
    ignored, and for a `DispersalKernel`, all other keywords are ignored.
- `formulation`: kernel formulation object holding the exact form of the kernal.
    Default [`ExponentialKernel`](@ref).
- `cellsize`: the cell size of the discretised kernal (i.e. simulation grid size).
    Default is 1.0.
- `distancemethod`: [`DistanceMethod`](@ref) object for calculating distance between cells.
    The default is [`CentroidToCentroid`](@ref).
"""
struct DispersalKernel{R,N,L,T,H<:Stencil{R,N,L,T},K,F,C,D<:DistanceMethod} <: Stencils.AbstractKernelStencil{R,N,L,T,H}
    stencil::H
    kernel::K
    formulation::F
    cellsize::C
    distancemethod::D
    function DispersalKernel(
        stencil::H, kernel, formulation::F, cellsize::C, distancemethod::D
    ) where {H<:Stencil{R,N,L,T},F,C,D<:DistanceMethod} where {R,N,L,T}
        if stencil isa DispersalKernel
            stencil
        else
            # Build the kernel matrix
            newkernel = scale(buildkernel(stencil, formulation, distancemethod, cellsize))
            new{R,N,L,T,H,typeof(newkernel),F,C,D}(
                stencil, newkernel, formulation, cellsize, distancemethod
            )
        end
    end
    function DispersalKernel{R,N,L,T,H,K,F,C,D}(
        stencil::H, kernel::K, formulation::F, cellsize::C, distancemethod::D
    ) where {R,N,L,T,H,K,F,C,D}
        new{R,N,L,T,H,K,F,C,D}(stencil, kernel, formulation, cellsize, distancemethod)
    end
end
function DispersalKernel(;
    radius=1,
    stencil=Window{radius}(),
    formulation=ExponentialKernel(),
    cellsize=1.0,
    distancemethod=CentroidToCentroid(),
)
    DispersalKernel(stencil, nothing, formulation, cellsize, distancemethod)
end
DispersalKernel{R}(; radius=R, kw...) where R = DispersalKernel(; radius=radius, kw...)

ConstructionBase.constructorof(::Type{<:DispersalKernel}) = DispersalKernel

function DG.Stencils.rebuild(n::DispersalKernel{R,N,L,<:Any,<:Any,K,F,C,D}, buffer) where {R,N,L,K,F,C,D}
    newstencil = Stencils.rebuild(stencil(n), buffer)
    DispersalKernel{R,N,L,eltype(newstencil),typeof(newstencil),K,F,C,D}(
        newstencil, kernel(n), formulation(n), cellsize(n), distancemethod(n)
    )
end

cellsize(stencil::DispersalKernel) = stencil.cellsize
distancemethod(stencil::DispersalKernel) = stencil.distancemethod
formulation(stencil::DispersalKernel) = stencil.formulation

# function buildkernel(window::Window{R}, formulation, distancemethod, cellsize) where R
#     # The radius doesn't include the center cell, so add it
#     S = 2R + 1
#     kernel = zeros(typeof(cellsize), S, S)
#     r1 = R + one(R)
#     # Paper: l. 97
#     for x = 0:R, y = 0:R
#         # Calculate the distance effect from the center cell to this cell
#         prob = dispersalprob(formulation, distancemethod, x, y, cellsize)
#         # Update the kernel value based on the formulation and distance
#         kernel[ x + r1,  y + r1] = prob
#         kernel[ x + r1, -y + r1] = prob
#         kernel[-x + r1,  y + r1] = prob
#         kernel[-x + r1, -y + r1] = prob
#     end
#     SMatrix{S,S}(kernel)
# end
function buildkernel(stencil::Stencil{<:Any,<:Any,L}, f, dm, cellsize) where L
    SVector{L}(Tuple(dispersalprob(f, dm, x, y, cellsize) for (x, y) in offsets(stencil)))
end

scale(x) = x ./ sum(x)
