"""
    DispersalKernel{R}(; 
        neighborhood=Window{R=1}(),
        formulation=ExponentialKernel(1.0), 
        cellsize=1.0, 
        distancemethod=CentroidToCentroid()
    )

A discretised two-dimensional dispersal kernel. May hold any `Neighborhood` object: the 
kernel will be built to match the shape, following the `formulation`, `cellsize`
and `distancemethod`.

# Keyword Arguments

- `neighborhood`: `Neighborhood` object specifying the range from the origin of the 
discretised dispersal kernal.
- `formulation`: kernel formulation object holding the exact form of the kernal.
- `cellsize`: the cell size of the discretised kernal (i.e. simulation grid size)
- `distancemethod`: method object for calculating distance between cells. The default is [`ExponentialKernel`](@ref).
    
`R` will be taken from the neighborhood, unless no neighborhood is provided.
"""
struct DispersalKernel{R,L,N<:Neighborhood{R,L},K,F,C,D} <: AbstractKernel{R,L}
    neighborhood::N
    kernel::K
    formulation::F
    cellsize::C
    distancemethod::D
end
function DispersalKernel(hood::N, kernel, formulation::F, cellsize::C, distancemethod::D
) where {N<:Neighborhood{R,L},F,C,D} where {R,L}
    # Build the kernel matrix
    newkernel = scale(buildkernel(hood, formulation, distancemethod, cellsize))
    DispersalKernel{R,L,N,typeof(newkernel),F,C,D}(
        hood, newkernel, formulation, cellsize, distancemethod
    )
end
function DispersalKernel(;
    neighborhood=Window{1}(),
    formulation=ExponentialKernel(),
    cellsize=1.0,
    distancemethod=CentroidToCentroid(),
)
    DispersalKernel(neighborhood, nothing, formulation, cellsize, distancemethod)
end
function DispersalKernel{R}(; neighborhood=Window{R}(), kw...) where R 
    DispersalKernel(; neighborhood=neighborhood, kw...)
end

ConstructionBase.constructorof(::Type{<:DispersalKernel}) = DispersalKernel 

function DG._setbuffer(n::DispersalKernel{R,L,<:Any,K,F,C,D}, buffer) where {R,L,K,F,C,D}
    newhood = DG._setbuffer(neighborhood(n), buffer)
    DispersalKernel{R,L,typeof(newhood),K,F,C,D}(
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
function buildkernel(window::Neighborhood{<:Any,L}, f, dm, cellsize) where L
    SVector{L}(Tuple(dispersalprob(f, dm, x, y, cellsize) for (x, y) in offsets(window)))
end

scale(x) = x ./ sum(x)

"""
Abstract supertype for methods of calculating distances and discretised dispersal 
probabilities between cells in a grid.
"""
abstract type DistanceMethod end

subsampling(method::DistanceMethod) = method.subsampling

"""
    CentroidToCentroid()

Calculates the discrete probability of dispersal between source and destination cell 
centroids. This is the naive method, but it will not handle low grid resolution well 
due to severe truncation.
"""
struct CentroidToCentroid <: DistanceMethod end

dispersalprob(f, ::CentroidToCentroid, x, y, cellsize) = f(sqrt(x^2 + y^2) * cellsize)

# """
#     CentroidToCentroid(subsampling)
#     CentroidToCentroid(; subsampling=10.0)
#
# Calculates probability of dispersal between source cell centroid and destination cell area.
# """
# @columns struct CentroidToArea <: DistanceMethod
    # Field        | Default | Flat | Bounds      | Description
    # subsampling::Int | 10.0    | true | (2.0, 40.0) | "Subsampling for brute-force integration"
# end

# dispersalprob(f, dm::CentroidToArea, x, y, cellsize) = error("not implemented yet")

"""
    AreaToCentroid(subsampling)
    AreaToCentroid(; subsampling=10.0)
    
Calculates the discrete probability of dispersal between source cell area and destination
centroid.
"""
Base.@kwdef struct AreaToCentroid{SS<:Number} <: DistanceMethod
    "Subsampling for numerical integration"
    subsampling::SS = Param(10.0; bounds=(2.0, 40.0))
end

@inline function dispersalprob(f, dm::AreaToCentroid, x, y, cellsize)
    prob = zero(cellsize)
    centerfirst = 1 / subsampling(dm) / 2 - 0.5
    centerlast = centerfirst * -1
    range = LinRange(centerfirst, centerlast, round(Int, subsampling(dm)))
    for j in range, i in range
        prob += sqrt((x + i)^2 + (y + j)^2) * cellsize |> f
    end
    prob / subsampling(dm)^2
end

"""
    AreaToArea(subsampling)
    AreaToArea(; subsampling=10.0)


Calculates the discrete probability of dispersal between source and destination based on 
cell areas.
"""
Base.@kwdef struct AreaToArea{SS<:Number} <: DistanceMethod
    subsampling::SS = Param(10.0; bounds=(2.0, 40.0))
end

@inline function dispersalprob(f, dm::AreaToArea, x, y, cellsize)
    prob = zero(cellsize)
    # Get the center point of the first cell (for both dimensions)
    centerfirst = 1 / subsampling(dm) / 2 - 0.5
    centerlast = centerfirst * -1
    range = LinRange(centerfirst, centerlast, round(Int, subsampling(dm)))
    for i in range, j in range
        for a in range, b in range
            prob += sqrt((x + i + a)^2 + (y + j + b)^2) * cellsize |> f
        end
    end
    prob / subsampling(dm)^4
end


"""
Abstract supertype for functors that calculate the probability density based on distance.

Concrete implementations must define functor methods with the form:

```julia
(k::SomeKernel)(distance) = # do something with `distance` and `k`"
```

Using an anonymous function would not allow rebuildable model parameters.
"""
abstract type KernelFormulation end

"""
    ExponentialKernel(λ)

Probability density function of distance ``d``.
```math
y = e^{-d/λ}
```
where λ is a shape parameter.
"""
@Base.kwdef struct ExponentialKernel{P} <: KernelFormulation
    "Parameter for adjusting spread of dispersal propability"
    λ::P = Param(1.0, bounds=(0.0, 2.0))
end
(f::ExponentialKernel)(distance) = exp(-distance / f.λ)
