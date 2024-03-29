"""
    DispersalKernel <: AbstractKernelNeighborhood

    DispersalKernel(; kw...)

Dispersal kernel for taking the dot product of the neighborhood and a matching
kernel of weights. May hold any `Neighborhood` object: the kernel
will be built to match the shape, using the `folumation`, `cellsize`
and `distancemethod`.

# Keyword Arguments

- `neighborhood`: Any DynamicGrids.jl `Neighborhood`, or an
    already constructed [`DispersalKernel`](@ref). Using this keyword means `radius` is
    ignored, and for a `DispersalKernel`, all other keywords are ignored.
- `neighborhood`: `Neighborhood` object specifying the range from the origin of the
    discretised dispersal kernal. Defaults to `Window(radius)`.
- `formulation`: kernel formulation object holding the exact form of the kernal.
    Default [`ExponentialKernel`](@ref).
- `cellsize`: the cell size of the discretised kernal (i.e. simulation grid size).
    Default is 1.0.
- `distancemethod`: [`DistanceMethod`](@ref) object for calculating distance between cells.
    The default is [`CentroidToCentroid`](@ref).
"""
struct DispersalKernel{R,N,L,H<:Neighborhood{R,N,L},K,F,C,D<:DistanceMethod} <: AbstractKernelNeighborhood{R,N,L,H}
    neighborhood::H
    kernel::K
    formulation::F
    cellsize::C
    distancemethod::D
    function DispersalKernel(
        hood::H, kernel, formulation::F, cellsize::C, distancemethod::D
    ) where {H<:Neighborhood{R,N,L},F,C,D<:DistanceMethod} where {R,N,L}
        if hood isa DispersalKernel
            hood
        else
            # Build the kernel matrix
            newkernel = scale(buildkernel(hood, formulation, distancemethod, cellsize))
            new{R,N,L,H,typeof(newkernel),F,C,D}(
                hood, newkernel, formulation, cellsize, distancemethod
            )
        end
    end
    function DispersalKernel{R,N,L,H,K,F,C,D}(
        hood::H, kernel::K, formulation::F, cellsize::C, distancemethod::D
    ) where {R,N,L,H,K,F,C,D}
        new{R,N,L,H,K,F,C,D}(hood, kernel, formulation, cellsize, distancemethod)
    end
end
function DispersalKernel(;
    radius=1,
    neighborhood=Window{radius}(),
    formulation=ExponentialKernel(),
    cellsize=1.0,
    distancemethod=CentroidToCentroid(),
)
    DispersalKernel(neighborhood, nothing, formulation, cellsize, distancemethod)
end
DispersalKernel{R}(; radius=R, kw...) where R = DispersalKernel(; radius=radius, kw...)

ConstructionBase.constructorof(::Type{<:DispersalKernel}) = DispersalKernel

function DG.setwindow(n::DispersalKernel{R,N,L,<:Any,K,F,C,D}, buffer) where {R,N,L,K,F,C,D}
    newhood = DG.setwindow(neighborhood(n), buffer)
    DispersalKernel{R,N,L,typeof(newhood),K,F,C,D}(
        newhood, kernel(n), formulation(n), cellsize(n), distancemethod(n)
    )
end

cellsize(hood::DispersalKernel) = hood.cellsize
distancemethod(hood::DispersalKernel) = hood.distancemethod
formulation(hood::DispersalKernel) = hood.formulation

function buildkernel(window::Window{R}, formulation, distancemethod, cellsize) where R
    # The radius doesn't include the center cell, so add it
    S = 2R + 1
    kernel = zeros(typeof(cellsize), S, S)
    r1 = R + one(R)
    # Paper: l. 97
    for x = 0:R, y = 0:R
        # Calculate the distance effect from the center cell to this cell
        prob = dispersalprob(formulation, distancemethod, x, y, cellsize)
        # Update the kernel value based on the formulation and distance
        kernel[ x + r1,  y + r1] = prob
        kernel[ x + r1, -y + r1] = prob
        kernel[-x + r1,  y + r1] = prob
        kernel[-x + r1, -y + r1] = prob
    end
    SMatrix{S,S}(kernel)
end
function buildkernel(window::Neighborhood{<:Any,<:Any,L}, f, dm, cellsize) where L
    SVector{L}(Tuple(dispersalprob(f, dm, x, y, cellsize) for (x, y) in offsets(window)))
end

scale(x) = x ./ sum(x)
