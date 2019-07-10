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

@columns struct DispersalKernel{F,K,C,I} <: AbstractDispersalKernel
    formulation::F | ExponentialKernel(1.0) | true  | _           | "Kernel formulation object"
    kernel::K      | []                     | false | _           | "Kernal matrix"
    cellsize::C    | 1.0                    | false | (0.0, 10.0) | "Simulation cell size"
    radius::I      | 3                      | false | (1, 10)     | "Kernel radius"

    DispersalKernel(formulation::F, init_kernel::K, cellsize::C, radius::I) where {F,K,C,I} = begin
        # Convert kenel the type of the init array
        inittype = typeof(init_kernel).name.wrapper
        kernel = inittype(build_dispersal_kernel(formulation, cellsize, radius))
        new{F,typeof(kernel),C,I}(formulation, kernel, cellsize, radius)
    end
end

build_dispersal_kernel(formulation, cellsize, r) = begin
    # The radius doesn't include the center cell, so add it
    sze = 2r + 1
    kernel = zeros(Float64, sze, sze)
    # Paper: l. 97
    for x = -r:r, y = -r:r
        # Calculate the distance from the center cell to this cell
        distance = sqrt(x^2 + y^2) * cellsize
        # Update the kernel value based on the formulation and distance
        kernel[x + r + one(x), y + r + one(y)] = 
            dispersalatdistance(formulation, distance)
    end
    # Normalise
    kernel ./= sum(kernel)
end

CellularAutomataBase.radius(hood::DispersalKernel) = hood.radius



abstract type AbstractKernelFormulation end

@description @limits @flattenable struct ExponentialKernel{P} <: AbstractKernelFormulation
    param::P    | true  | (0.0, 100.0) | "Parameter for adjusting kernel spread"
end

# Paper: l. 96 
dispersalatdistance(f::ExponentialKernel, distance) = exp(-distance / f.param)

