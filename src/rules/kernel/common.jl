"Extends AbstractNeighborhood for for dispersal kernel neighborhoods"
abstract type AbstractDispersalKernel{R} <: AbstractRadialNeighborhood{R} end

@mix @columns struct Kernel{N}
    neighborhood::N | DispersalKernel{3}() | true  | _ | "Normalised proportions of dispersal to surrounding cells"
end

"""
A neighborhood built from a dispersal kernel function and a cell radius.
Can be built directly by passing in the array, radius and overflow
arguments, but preferably use the keyword constructor to build the array from
a dispersal kernel function.
$(FIELDDOCTABLE)
"""
@columns struct DispersalKernel{R,F,K,C,D} <: AbstractDispersalKernel{R}
    formulation::F    | ExponentialKernel(1.0) | true  | _           | "Kernel formulation object"
    kernel::K         | nothing                | false | _           | "Kernal matrix"
    cellsize::C       | 1.0                    | false | (0.0, 10.0) | "Simulation cell size"
    distancemethod::D | CentroidToCentroid()   | false | _           | "Method for calculating distance between cells"

    DispersalKernel{R}(; kwargs...) where {R,F,K,C,D} = begin
        args = FieldDefaults.insert_kwargs(kwargs, DispersalKernel)
        DispersalKernel{R,typeof.(args)...}(args...)
    end
    DispersalKernel{R}(formulation::F, kernel::K, cellsize::C, distancemethod::D) where {R,F,K,C,D} =
        DispersalKernel{R,F,K,C,D}(formulation, kernel, cellsize, distancemethod)
    DispersalKernel{R,F,K,C,D}(formulation::F, kernel::K, cellsize::C, distancemethod::D) where {R,F,K,C,D} = begin
        # Convert kenel the type of the init array
        kernel = scale(buildkernel(formulation, distancemethod, cellsize, R))
        kernel = K <: Nothing ? kernel : K(kernel)
        new{R,F,typeof(kernel),C,D}(formulation, kernel, cellsize, distancemethod)
    end
end

buildkernel(formulation, distancemethod, cellsize, r) = begin
    # The radius doesn't include the center cell, so add it
    sze = 2r + 1
    kernel = zeros(typeof(cellsize), sze, sze)
    r1 = r + one(r)
    # Paper: l. 97
    for x = 0:r, y = 0:r
        # Calculate the distance from the center cell to this cell
        prob = dispersalprob(formulation, distancemethod, x, y, cellsize)
        # Update the kernel value based on the formulation and distance
        kernel[x + r1, y + r1] = prob
        kernel[x + r1, -y + r1] = prob
        kernel[-x + r1, y + r1] = prob
        kernel[-x + r1, -y + r1] = prob
    end
    kernel
end

scale(x) = x ./= sum(x)


# build_dispersal_kernel(formulation, distancemethod, cellsize, r) = begin
#     # The radius doesn't include the center cell, so add it
#     sze = 2r + 1
#     kernels = [zeros(Float64, sze, sze) for thread in 1:Threads.nthreads()]
#     r1 = r + one(r)

#     # Calculate one quadrant
#     Threads.@threads for x = 0:r
#         kernel = kernels[Threads.threadid()]
#         for y = 0:r
#             kernel[x + r1, y + r1] = dispersalprob(formulation, distancemethod, x, y, cellsize)
#         end
#     end
#     kernel = sum(kernels)
#     println()
#     println(length(kernels))
#     display(kernels)
#     println()
#     # Update the other 3 quadrants
#     for x = 0:r, y = 0:r
#         # Just sum the threads, most are zero
#         prob = kernel[x + r1, y + r1]
#         kernel[x + r1, -y + r1] = prob
#         kernel[-x + r1, y + r1] = prob
#         kernel[-x + r1, -y + r1] = prob
#     end
#     # Normalise
#     kernel ./= sum(kernel)
# end

DynamicGrids.radius(hood::DispersalKernel{R}) where R = R

Flatten.constructor_of(::Type{<:DispersalKernel{R}}) where R = DispersalKernel{R}


"""
Methods for calculating distances and dispersal probabilities between cells in a grid.
"""
abstract type AbstractDistanceMethod end

"""
Calculates probability of dispersal between source and destination cell centroids
This is the obvious, naive method, but it will not handle low grid resolution well.
"""
struct CentroidToCentroid <: AbstractDistanceMethod end

dispersalprob(f, ::CentroidToCentroid, x, y, cellsize) = sqrt(x^2 + y^2) * cellsize |> f

"""
Calculates probability of dispersal between source cell centroid and destination cell area.
"""
@columns struct CentroidToArea <: AbstractDistanceMethod
    subsample::Int | 10.0 | true | (2.0, 40.0) | "Subsampling for brute-force integration"
end

dispersalprob(f, dm::CentroidToArea, x, y, cellsize) = error("not implemented yet")

"""
Calculates probability of dispersal between source cell area and destination centroid.
"""
@columns struct AreaToCentroid <: AbstractDistanceMethod
    subsample::Int | 10.0 | true | (2.0, 40.0) | "Subsampling for brute-force integration"
end
AreaToCentroid(subsample::Float64) = AreaToCentroid(round(Int, subsample))

@inline dispersalprob(f, dm::AreaToCentroid, x, y, cellsize) = begin
    prob = zero(cellsize)
    centerfirst = 1 / dm.subsample / 2 - 0.5
    centerlast = centerfirst * -1
    range = LinRange(centerfirst, centerlast, dm.subsample)
    for j in range, i in range
        prob += sqrt((x + i)^2 + (y + j)^2) * cellsize |> f
    end
    prob / dm.subsample^2
end

"""
Calculates probability of dispersal between source and destination cell areas.
"""
struct AreaToArea <: AbstractDistanceMethod
    subsample::Int
end
AreaToArea(subsample::Float64) = AreaToArea(round(Int, subsample))

@inline dispersalprob(f, dm::AreaToArea, x, y, cellsize) = begin
    prob = zero(cellsize)
    # Get the center point of the first cell (for both dimensions)
    centerfirst = 1 / dm.subsample / 2 - 0.5
    centerlast = centerfirst * -1
    range = LinRange(centerfirst, centerlast, dm.subsample)
    for i in range, j in range
        for a in range, b in range
            prob += sqrt((x + i + a)^2 + (y + j + b)^2) * cellsize |> f
        end
    end
    prob / dm.subsample^4
end


"""
Functors for calculating the probability of dispersal between two points.
"""
abstract type AbstractKernelFormulation end

"""
Probability of dispersal with a negatitve exponential relationship to distance.
"""
@description @limits @flattenable struct ExponentialKernel{P} <: AbstractKernelFormulation
    λ::P    | true  | (0.0, 2.0) | "Parameter for adjusting spread of dispersal propability"
end
(f::ExponentialKernel)(distance) = exp(-distance / f.λ)
