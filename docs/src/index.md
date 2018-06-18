# Dispersal.jl

```@meta
CurrentModule = Dispersal
```

```@docs
Dispersal
```

## Example

This is the general pattern of running model in Dispersal.jl. 

We declare a suitability layer and initial conditions. Then we run a dispersal
simulation that combines local and jump dispersal, for three timesteps.

```@example
using Cellular
using Dispersal

# add a 0.0-1.0 scaled raster array that represents habitat suitability
suitability = [0.5 0.0 0.3 0.0 0.0 0.3 0.0;
               0.0 0.2 0.8 0.0 0.9 0.6 0.0;
               0.0 0.5 1.0 1.0 1.0 0.0 0.0;
               0.0 0.0 1.0 1.0 0.7 0.0 0.0;
               0.0 0.0 0.6 1.0 0.0 0.0 0.0;
               0.0 0.1 0.6 0.0 0.0 0.0 0.0;
               0.0 0.0 0.0 0.0 0.0 0.0 0.7]

# Make an init array the same size as your suitability layer, and seed it
init = [0  0  0  0  0  0  0;
        0  0  0  0  0  0  0;
        0  0  0  0  0  0  0;
        0  0  0  0  0  0  0;
        0  0  1  0  0  0  0;
        0  0  0  0  0  0  0;
        0  0  0  0  0  0  0]

# You also could do this with:
# init = zeros(suitability); init[5, 3] = 1


# Define a dispersal kernel function

f = d -> e^-d

# Define the neighborhood, using the dispersal kernel and a radius
hood = DispersalNeighborhood(; f=f, radius=2)

# Define additional raster layers
layers = SuitabilityLayer(suitability)

# Define disersal modules
localdisp = LocalDispersal(layers=layers, neighborhood=hood)
jumpdisp = JumpDispersal(layers=layers)

# Set the output type
output = ArrayOutput(init)

# Run the simulation
sim!(output, (localdisp, jumpdisp), init; time=1:3) 

output.frames[3]
```

## Models

```@docs
AbstractLocalDispersal
LocalDispersal
AbstractJumpDispersal
JumpDispersal
AbstractHumanDispersal
HumanDispersal
```

These modules all trigger methods of `rule`:

```@docs
rule
```


## Neighborhoods

Extend Cellular.AbstractNeighborhood

```@docs
AbstractDispersalNeighborhood
DispersalNeighborhood
```

### Customisation

```@docs
neighbors
build_dispersal_kernel
pressure
```

## Layers

Layers are a concept not present in Cellular.jl. They provide 
overlay grids of additional information about dispersal potential.
Raw arrays are wrapped in a type that determines how they will be used.

```@docs
AbstractLayer 
AbstractSuitabilityLayer 
SuitabilityLayer 
AbstractSuitabilitySequence 
SuitabilitySequence
HumanLayer
AbstractLayers
Layers
```

### Customisation

```@docs
suitability
seq_interpolate
get_cell
```
