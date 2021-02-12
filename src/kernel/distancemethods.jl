"""
    DistanceMethod

Abstract supertype for methods of calculating distances and discretised dispersal 
probabilities between cells in a grid.
between cells in a grid.

Distance calculation methods include:
- [`CentroidToCentroid`](@ref)
- [`AreaToArea`](@ref)
- [`AreaToCentroid`](@ref)

Which are adapted from: "On the approximation of continuous dispersal kernels 
in discrete-space models, Joseph D. Chipperfield et al 2011"

The `CentroidToArea` method has not been implemented.
"""
abstract type DistanceMethod end

subsampling(method::DistanceMethod) = method.subsampling

"""
    CentroidToCentroid <: DistanceMethod

    CentroidToCentroid()

Calculates the discrete probability of dispersal between source and destination cell 
centroids. This is the naive method, but it will not handle low grid resolution well 
due to severe truncation.

See: "On the approximation of continuous dispersal kernels in discrete-space models, 
Joseph D. Chipperfield et al 2011"
"""
struct CentroidToCentroid <: DistanceMethod end

dispersalprob(f, ::CentroidToCentroid, x, y, cellsize) = f(sqrt(x^2 + y^2) * cellsize)

"""
    AreaToCentroid <: DistanceMethod

    AreaToCentroid(subsampling)
    AreaToCentroid(; subsampling=10.0)
    
Calculates the discrete probability of dispersal between source cell area and destination
centroid.

See: "On the approximation of continuous dispersal kernels in discrete-space models, 
Joseph D. Chipperfield et al 2011"
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
    AreaToArea <: DistanceMethod

    AreaToArea(subsampling)
    AreaToArea(; subsampling=10.0)

Calculates the discrete probability of dispersal between source and destination based on 
cell areas.

See: "On the approximation of continuous dispersal kernels in discrete-space models, 
Joseph D. Chipperfield et al 2011"
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
