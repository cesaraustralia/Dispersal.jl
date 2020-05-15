@mix @columns struct Kernel{NH}
    # Field          | Default              | Flat | Bounds | Description
    neighborhood::NH | DispersalKernel{3}() | true | _      | "Normalised proportions of dispersal to surrounding cells"
end

"""
    DispersalKernel{Radius}(formulation::F, kernel::K, cellsize::C, distancemethod::D)
    DispersalKernel{Radius}(; formulation=ExponentialKernel(1.0), cellsize=1.0, distancemethod=CentroidToCentroid())

arguments, but preferably use the keyword constructor to build the array from
a dispersal kernel function.
$(FIELDDOCTABLE)
"""
@columns struct DispersalKernel{R,F,K,C,D} <: AbstractRadialNeighborhood{R}
    # Field           | Default                | Flat  | Bounds      | Description
    formulation::F    | ExponentialKernel(1.0) | true  | _           | "Kernel formulation object"
    kernel::K         | nothing                | false | _           | "Kernal matrix"
    cellsize::C       | 1.0                    | false | (0.0, 10.0) | "Simulation cell size"
    distancemethod::D | CentroidToCentroid()   | false | _           | "Method for calculating distance between cells"

    DispersalKernel{R,F,K,C,D}(formulation::F, kernel::K, cellsize::C, distancemethod::D) where {R,F,K,C,D} = begin
        # Convert kenel the type of the init array
        kernel = scale(buildkernel(formulation, distancemethod, cellsize, R))
        kernel = K <: Nothing ? kernel : K(kernel)
        new{R,F,typeof(kernel),C,D}(formulation, kernel, cellsize, distancemethod)
    end
end
DispersalKernel{R}(; kwargs...) where {R,F,K,C,D} = begin
    args = FieldDefaults.insert_kwargs(kwargs, DispersalKernel)
    DispersalKernel{R,typeof.(args)...}(args...)
end
DispersalKernel{R}(formulation::F, kernel::K, cellsize::C, distancemethod::D) where {R,F,K,C,D} =
    DispersalKernel{R,F,K,C,D}(formulation, kernel, cellsize, distancemethod)

ConstructionBase.constructorof(::Type{<:DispersalKernel{R}}) where R = DispersalKernel{R}

DynamicGrids.radius(hood::DispersalKernel{R}) where R = R
kernel(hood::DispersalKernel) = hood.kernel
formulation(hood::DispersalKernel) = hood.formulation

@inline applykernel(hood::DispersalKernel, buf) =
    @inbounds return buf ⋅ kernel(hood)


buildkernel(formulation, distancemethod, cellsize, r) = begin
    # The radius doesn't include the center cell, so add it
    sze = 2r + 1
    kernel = zeros(typeof(cellsize), sze, sze)
    r1 = r + one(r)
    # Paper: l. 97
    for x = 0:r, y = 0:r
        # Calculate the distance effect from the center cell to this cell
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


"""
Abstract supertype for methods of calculating distances and dispersal probabilities
between cells in a grid.
"""
abstract type DistanceMethod end

subsampling(method::DistanceMethod) = method.subsampling

"""
    CentroidToCentroid()

Calculates probability of dispersal between source and destination cell centroids
This is the obvious, naive method, but it will not handle low grid resolution well.
"""
struct CentroidToCentroid <: DistanceMethod end

dispersalprob(f, ::CentroidToCentroid, x, y, cellsize) = 
    sqrt(x^2 + y^2) * cellsize |> f

# """
#     CentroidToCentroid(subsampling)
#     CentroidToCentroid(; subsampling=10.0)
# 
# Calculates probability of dispersal between source cell centroid and destination cell area.
# 
# $(FIELDDOCTABLE)
# """
# @columns struct CentroidToArea <: DistanceMethod
    # Field        | Default | Flat | Bounds      | Description
    # subsampling::Int | 10.0    | true | (2.0, 40.0) | "Subsampling for brute-force integration"
# end

# dispersalprob(f, dm::CentroidToArea, x, y, cellsize) = error("not implemented yet")

"""
    AreaToCentroid(subsampling)
    AreaToCentroid(; subsampling=10.0)

Calculates probability of dispersal between source cell area and destination centroid.

$(FIELDDOCTABLE)
"""
@columns struct AreaToCentroid <: DistanceMethod
    # Field        | Default | Flat | Bounds      | Description
    subsampling::Int | 10.0    | true | (2.0, 40.0) | "Subsampling for brute-force integration"
end
AreaToCentroid(subsampling::AbstractFloat) = AreaToCentroid(round(Int, subsampling))

@inline dispersalprob(f, dm::AreaToCentroid, x, y, cellsize) = begin
    prob = zero(cellsize)
    centerfirst = 1 / subsampling(dm) / 2 - 0.5
    centerlast = centerfirst * -1
    range = LinRange(centerfirst, centerlast, subsampling(dm))
    for j in range, i in range
        prob += sqrt((x + i)^2 + (y + j)^2) * cellsize |> f
    end
    prob / subsampling(dm)^2
end

"""
    AreaToArea(subsampling)
    AreaToArea(; subsampling=10.0)


Calculates probability of dispersal between source and destination cell areas.

$(FIELDDOCTABLE)
"""
struct AreaToArea <: DistanceMethod
    subsampling::Int
end
AreaToArea(subsampling::AbstractFloat) = AreaToArea(round(Int, subsampling))

@inline dispersalprob(f, dm::AreaToArea, x, y, cellsize) = begin
    prob = zero(cellsize)
    # Get the center point of the first cell (for both dimensions)
    centerfirst = 1 / subsampling(dm) / 2 - 0.5
    centerlast = centerfirst * -1
    range = LinRange(centerfirst, centerlast, subsampling(dm))
    for i in range, j in range
        for a in range, b in range
            prob += sqrt((x + i + a)^2 + (y + j + b)^2) * cellsize |> f
        end
    end
    prob / subsampling(dm)^4
end


"""
Abstract supertype for functors for calculating the probability of 
dispersal between two points.

Concrete implementations must define functor methods with the form:
```julia
(k::SomeKernel)(x) = ? # do something with `x` and `k`"
```
"""
abstract type KernelFormulation end

"""
    ExponentialKernel(λ)

Probability of dispersal with a negatitve exponential relationship to distance.

$(FIELDDOCTABLE)
"""
@description @limits @flattenable struct ExponentialKernel{P} <: KernelFormulation
    λ::P    | true  | (0.0, 2.0) | "Parameter for adjusting spread of dispersal propability"
end
(f::ExponentialKernel)(distance) = exp(-distance / f.λ)
