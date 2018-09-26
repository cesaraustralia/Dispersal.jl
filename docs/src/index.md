# Dispersal.jl

```@meta
CurrentModule = Dispersal
```

```@docs
Dispersal
```

```@contents
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

# Or use ArchGDAL to import a raster imaage:
using ArchGDAL
function readtiff(file)
    img = ArchGDAL.registerdrivers() do
        ArchGDAL.read(file) do dataset
            ArchGDAL.read(dataset)
        end
    end
    # Scale values to a maximum of 1.0
    img ./= maximum(img)
    # Remove values below 0.0
    img .= max.(0.0, img)
    # Transpose: fix until ArchGDAL does this automatically (soon)
    img = img[:,:,1]'
end

suitability = readtiff("suitability.tif")

# You could also crop the raster, this function cuts out Australia from a world map:
cropaust(x) = x[950:1350, 3100:3600]

austsuitability = cropaust(suitability)

# Define a dispersal kernel function

f = d -> e^-d

# Define the neighborhood, using the dispersal kernel and a radius
hood = DispersalNeighborhood(; f=f, radius=2, init=init)

# Define raster layers. This tuple syntax lets you use layer arrays on GPUs but still
label their purpose, and dispatch on the layer type.
layers = (SuitabilityLayer(), suitability)

# Define disersal modules
localdisp = InwardsLocalDispersal(neighborhood=hood)
jumpdisp = JumpDispersal()

# Set the output type
output = ArrayOutput(init)

# Run the simulation
sim!(output, Models(localdisp, jumpdisp), init, layers; time=1:3) 

output.frames[3]
```

## Models

### Extending AbstractModel

```@docs
AbstractInwardsLocalDispersal
InwardsLocalDispersal
```

### Extending AbstractPartialModel

These models trigger rules that only partially update the grid.
Many dispersal functions only operate on cells that are currently occupied.

```@docs
AbstractOutwardsLocalDispersal
OutwardsLocalDispersal
AbstractJumpDispersal
JumpDispersal
AbstractHumanDispersal
HumanDispersal
```

Models trigger methods of [`rule`](@ref). 
Custom models you create will in most cases need a custom `rule` method.

```@docs
rule
```


## Neighborhoods

Extend Cellular.AbstractNeighborhood, and add `neighbors()` methods.

### Types and Constructors

```@docs
AbstractDispersalNeighborhood
DispersalNeighborhood
DispersalNeighborhood(; dir=:inwards, f=exponential, param=1.0, init=[], cellsize=1.0, radius=3, overflow=Skip())
```

### Methods

```@docs
neighbors
neighbors(hood::DispersalNeighborhood, model, state, row, col, t, source, dest, args...) = begin
pressure
```

## Layers

Layers are a concept not present in Cellular.jl. They provide 
overlay grids of additional information about dispersal potential.

Like models, than can be combined arbitrarily, in this case in a tuple. Methods
loop through all layers to return a scalar that is the product of their
outputs. Unrelated layers are ignored.

### Types

```@docs
AbstractLayer 
AbstractSuitabilityLayer 
SuitabilityLayer 
AbstractSuitabilitySequence 
SuitabilitySequence
HumanLayer
```

### Methods 

```@docs
suitability
human_impact
sequence_interpolate
cyclic
```
