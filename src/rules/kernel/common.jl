"Extends AbstractNeighborhood for for dispersal kernel neighborhoods"
abstract type AbstractDispersalKernel{R} <: AbstractNeighborhood{R} end

@mix @columns struct Kernel{N}
    neighborhood::N | DispersalKernel{3}() | true  | _ | "Neighborhood to disperse to or from"
end

"""
A neighborhood built from a dispersal kernel function and a cell radius.
Can be built directly by passing in the array, radius and overflow
arguments, but preferably use the keyword constructor to build the array from
a dispersal kernel function.
$(FIELDDOCTABLE)
"""

@columns struct DispersalKernel{R,F,K,C} <: AbstractDispersalKernel{R}
    formulation::F | ExponentialKernel(1.0) | true  | _           | "Kernel formulation object"
    kernel::K      | nothing                | false | _           | "Kernal matrix"
    cellsize::C    | 1.0                    | false | (0.0, 10.0) | "Simulation cell size"

    DispersalKernel{R}(; kwargs...) where {R,F,K,C} = begin
        kwargs = FieldDefaults.insert_kwargs(kwargs, DispersalKernel)
        DispersalKernel{R,typeof.(kwargs)...}(kwargs...)
    end
    DispersalKernel{R}(formulation::F, kernel::K, cellsize::C) where {R,F,K,C} =
        DispersalKernel{R,F,K,C}(formulation, kernel, cellsize)
    DispersalKernel{R,F,K,C}(formulation::F, kernel::K, cellsize::C) where {R,F,K,C} = begin
        # Convert kenel the type of the init array
        kernel = build_dispersal_kernel(formulation, CentroidToCentroid(), cellsize, R)
        kernel = K <: Nothing ? kernel : K(kernel)
        new{R,F,typeof(kernel),C}(formulation, kernel, cellsize)
    end
end


build_dispersal_kernel(formulation, distancemethod, cellsize, r) = begin
    # The radius doesn't include the center cell, so add it
    sze = 2r + 1
    kernel = zeros(Float64, sze, sze)
    # Paper: l. 97
    for x = -r:r, y = -r:r
        # Calculate the distance from the center cell to this cell
        dist = distance(distancemethod, x, y, cellsize)
        # Update the kernel value based on the formulation and distance
        kernel[x + r + one(x), y + r + one(y)] =
            dispersalatdistance(formulation, dist)
    end
    # Normalise
    kernel ./= sum(kernel)
end

CellularAutomataBase.radius(hood::DispersalKernel{R}) where R = R

Flatten.constructor_of(::Type{<:DispersalKernel{R}}) where R = DispersalKernel{R}


abstract type AbstractDistanceMethod end


struct CentroidToCentroid <: AbstractDistanceMethod end

distance(::CentroidToCentroid, x, y, cellsize) = sqrt(x^2 + y^2) * cellsize

struct CentroidToArea <: AbstractDistanceMethod end

distance(::CentroidToArea, x, y, cellsize) = error("not implemented yet")

struct AreaToArea <: AbstractDistanceMethod end

distance(::AreaToArea, x, y, cellsize) = error("not implemented yet")

struct AreaToCentroid <: AbstractDistanceMethod end

distance(::AreaToCentroid, x, y, cellsize) = error("not implemented yet")



abstract type AbstractKernelFormulation end

@description @limits @flattenable struct ExponentialKernel{P} <: AbstractKernelFormulation
    λ::P    | true  | (0.0, 3.0) | "Parameter for adjusting spread of dispersal propability"
end

# Paper: l. 96
dispersalatdistance(f::ExponentialKernel, distance) = exp(-distance / f.λ)
