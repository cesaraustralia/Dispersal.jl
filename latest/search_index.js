var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Dispersal",
    "title": "Dispersal",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Dispersal",
    "page": "Dispersal",
    "title": "Dispersal",
    "category": "module",
    "text": "This package extends Cellular.jl\n\nDispersal.jl provides a range of dispersal modules that can be combined to build grid-based organism dispersal simulations.\n\nThe design provides a solid framework while allowing customisation of any  aspect. A model may start with the defaults and formulations provided,  but incrementally customise them for a particular use-case, to any level of detail. \n\nAdditionally, modules, outputs, neighborhoods provided by Cellular.jl or other  packages that extend it may be incorporated into a simulaiton.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.jl-1",
    "page": "Dispersal",
    "title": "Dispersal.jl",
    "category": "section",
    "text": "CurrentModule = DispersalDispersal"
},

{
    "location": "index.html#Example-1",
    "page": "Dispersal",
    "title": "Example",
    "category": "section",
    "text": "This is the general pattern of running model in Dispersal.jl. We declare a suitability layer and initial conditions. Then we run a dispersal simulation that combines local and jump dispersal, for three timesteps.using Cellular\nusing Dispersal\n\n# add a 0.0-1.0 scaled raster array that represents habitat suitability\nsuitability = [0.5 0.0 0.3 0.0 0.0 0.3 0.0;\n               0.0 0.2 0.8 0.0 0.9 0.6 0.0;\n               0.0 0.5 1.0 1.0 1.0 0.0 0.0;\n               0.0 0.0 1.0 1.0 0.7 0.0 0.0;\n               0.0 0.0 0.6 1.0 0.0 0.0 0.0;\n               0.0 0.1 0.6 0.0 0.0 0.0 0.0;\n               0.0 0.0 0.0 0.0 0.0 0.0 0.7]\n\n# Make an init array the same size as your suitability layer, and seed it\ninit = [0  0  0  0  0  0  0;\n        0  0  0  0  0  0  0;\n        0  0  0  0  0  0  0;\n        0  0  0  0  0  0  0;\n        0  0  1  0  0  0  0;\n        0  0  0  0  0  0  0;\n        0  0  0  0  0  0  0]\n\n# You also could do this with:\n# init = zeros(suitability); init[5, 3] = 1\n\n\n# Define a dispersal kernel function\n\nf = d -> e^-d\n\n# Define the neighborhood, using the dispersal kernel and a radius\nhood = DispersalNeighborhood(; f=f, radius=2)\n\n# Define additional raster layers\nlayers = SuitabilityLayer(suitability)\n\n# Define disersal modules\nlocaldisp = InwardsLocalDispersal(layers=layers, neighborhood=hood)\njumpdisp = JumpDispersal(layers=layers)\n\n# Set the output type\noutput = ArrayOutput(init)\n\n# Run the simulation\nsim!(output, (localdisp, jumpdisp), init; time=1:3) \n\noutput.frames[3]"
},

{
    "location": "index.html#Models-1",
    "page": "Dispersal",
    "title": "Models",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Dispersal.AbstractInwardsLocalDispersal",
    "page": "Dispersal",
    "title": "Dispersal.AbstractInwardsLocalDispersal",
    "category": "type",
    "text": "abstract AbstractInwardsLocalDispersal <: Cellular.AbstractModel\n\n\n\n"
},

{
    "location": "index.html#Dispersal.InwardsLocalDispersal",
    "page": "Dispersal",
    "title": "Dispersal.InwardsLocalDispersal",
    "category": "type",
    "text": "Local dispersal within a DispersalNeighborhood or other neighborhoods. Inwards dispersal calculates dispersal to the current cell from cells in the neighborhood.\n\n\n\n"
},

{
    "location": "index.html#Extending-AbstractModel-1",
    "page": "Dispersal",
    "title": "Extending AbstractModel",
    "category": "section",
    "text": "AbstractInwardsLocalDispersal\nInwardsLocalDispersal"
},

{
    "location": "index.html#Dispersal.AbstractOutwardsLocalDispersal",
    "page": "Dispersal",
    "title": "Dispersal.AbstractOutwardsLocalDispersal",
    "category": "type",
    "text": "abstract AbstractOutwardsLocalDispersal <: Cellular.AbstractPartialModel\n\n\n\n"
},

{
    "location": "index.html#Dispersal.OutwardsLocalDispersal",
    "page": "Dispersal",
    "title": "Dispersal.OutwardsLocalDispersal",
    "category": "type",
    "text": "Local dispersal within a DispersalNeighborhood\n\nOutwards dispersal calculates dispersal from the current cell to cells  in its neighborhood. This should be more efficient than inwards  dispersal when a small number of cells are occupied, but less efficient when a large  proportion of the grid is occupied.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.AbstractJumpDispersal",
    "page": "Dispersal",
    "title": "Dispersal.AbstractJumpDispersal",
    "category": "type",
    "text": "abstract AbstractJumpDispersal <: Cellular.AbstractPartialModel\n\n\n\n"
},

{
    "location": "index.html#Dispersal.JumpDispersal",
    "page": "Dispersal",
    "title": "Dispersal.JumpDispersal",
    "category": "type",
    "text": "Jump dispersal within a DispersalNeighborhood] or other neighborhoods.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.AbstractHumanDispersal",
    "page": "Dispersal",
    "title": "Dispersal.AbstractHumanDispersal",
    "category": "type",
    "text": "abstract AbstractHumanDispersal <: Cellular.AbstractPartialModel\n\n\n\n"
},

{
    "location": "index.html#Dispersal.HumanDispersal",
    "page": "Dispersal",
    "title": "Dispersal.HumanDispersal",
    "category": "type",
    "text": "Human dispersal model.\n\n\n\n"
},

{
    "location": "index.html#Cellular.rule",
    "page": "Dispersal",
    "title": "Cellular.rule",
    "category": "function",
    "text": "rule(model::AbstractInwardsLocalDispersal, state, index, t, args...)\n\nRuns rule for of AbstractInwardsLocalDispersal dispersal. \n\nThe current cell is invaded if there is pressure from surrounding cells and  suitable habitat. Otherwise it keeps its current state.\n\n\n\nrule(model::AbstractOutwardsLocalDispersal, state, index, t, source, dest, args...)\n\nRuns rule for of AbstractOutwardsLocalDispersal dispersal. \n\nSurrounding cells are invaded if the current cell is occupied and they have  suitable habitat. Otherwise they keeps their current state.\n\n\n\nrule(model::AbstractJumpDispersal, state, index, t, source, dest, args...)\n\nLong range rule for AbstractJumpDispersal. A random cell within the spotrange is invaded if it is suitable.\n\n\n\nrule(model::AbstractHumanDispersal, state, index, t, source, dest, args...)\n\nSimulates human dispersal, weighting dispersal probability based on human population in the source cell.\n\n\n\n"
},

{
    "location": "index.html#Extending-AbstractPartialModel-1",
    "page": "Dispersal",
    "title": "Extending AbstractPartialModel",
    "category": "section",
    "text": "These models trigger rules that only partially update the grid. Many dispersal functions only operate on cells that are currently occupied.AbstractOutwardsLocalDispersal\nOutwardsLocalDispersal\nAbstractJumpDispersal\nJumpDispersal\nAbstractHumanDispersal\nHumanDispersalModels trigger methods of rule.  Custom models you create will in most cases need a custom rule method.rule"
},

{
    "location": "index.html#Neighborhoods-1",
    "page": "Dispersal",
    "title": "Neighborhoods",
    "category": "section",
    "text": "Extend Cellular.AbstractNeighborhood, and add neighbors() methods."
},

{
    "location": "index.html#Dispersal.AbstractDispersalNeighborhood",
    "page": "Dispersal",
    "title": "Dispersal.AbstractDispersalNeighborhood",
    "category": "type",
    "text": "abstract AbstractDispersalNeighborhood <: Cellular.AbstractNeighborhood\n\n\n\n"
},

{
    "location": "index.html#Dispersal.DispersalNeighborhood",
    "page": "Dispersal",
    "title": "Dispersal.DispersalNeighborhood",
    "category": "type",
    "text": "struct DispersalNeighborhood{K, S} <: Dispersal.AbstractDispersalNeighborhood\n\ndispkernel\nradius\noverflow\n\n\n\n"
},

{
    "location": "index.html#Dispersal.DispersalNeighborhood-Tuple{}",
    "page": "Dispersal",
    "title": "Dispersal.DispersalNeighborhood",
    "category": "method",
    "text": "DispersalNeighborhood(; f=d -> exponential(d, 1), radius=3, overflow=Skip())\n\nConstructor for neighborhoods, using a dispersal kernel function and a cell radius.\n\nKeyword Arguments:\n\nf::Function: any function that accepts a Number argument and returns a propbability between 0.0 and 1.0\nradius::Integer: a positive integer\n\n\n\n"
},

{
    "location": "index.html#Types-and-Constructors-1",
    "page": "Dispersal",
    "title": "Types and Constructors",
    "category": "section",
    "text": "AbstractDispersalNeighborhood\nDispersalNeighborhood\nDispersalNeighborhood(; f=d -> exponential(d, 1), radius=3, overflow=Skip())"
},

{
    "location": "index.html#Methods-1",
    "page": "Dispersal",
    "title": "Methods",
    "category": "section",
    "text": "neighbors\nneighbors(hood::DispersalNeighborhood, state, index, t, source, args...)\npressure"
},

{
    "location": "index.html#Layers-1",
    "page": "Dispersal",
    "title": "Layers",
    "category": "section",
    "text": "Layers are a concept not present in Cellular.jl. They provide  overlay grids of additional information about dispersal potential.Like models, than can be combined arbitrarily in tuples. Methods loop through all relevant layers to return a scalar that is the product of their outputs."
},

{
    "location": "index.html#Dispersal.AbstractLayer",
    "page": "Dispersal",
    "title": "Dispersal.AbstractLayer",
    "category": "type",
    "text": "abstract AbstractLayer\n\n\n\n"
},

{
    "location": "index.html#Dispersal.AbstractSuitabilityLayer",
    "page": "Dispersal",
    "title": "Dispersal.AbstractSuitabilityLayer",
    "category": "type",
    "text": "abstract AbstractSuitabilityLayer <: Dispersal.AbstractLayer\n\n\n\n"
},

{
    "location": "index.html#Dispersal.SuitabilityLayer",
    "page": "Dispersal",
    "title": "Dispersal.SuitabilityLayer",
    "category": "type",
    "text": "struct SuitabilityLayer{T} <: Dispersal.AbstractSuitabilityLayer\n\ndata\nAny 2-dimensional AbstractArray matching the coordinates of the init array\n\n\n\n"
},

{
    "location": "index.html#Dispersal.AbstractSuitabilitySequence",
    "page": "Dispersal",
    "title": "Dispersal.AbstractSuitabilitySequence",
    "category": "type",
    "text": "abstract AbstractSuitabilitySequence <: Dispersal.AbstractSuitabilityLayer\n\n\n\n"
},

{
    "location": "index.html#Dispersal.SuitabilitySequence",
    "page": "Dispersal",
    "title": "Dispersal.SuitabilitySequence",
    "category": "type",
    "text": "struct SuitabilitySequence{T, D} <: Dispersal.AbstractSuitabilitySequence\n\ntimespan\nThe timespan of each layer in the sequence.\ndata\nEither an Array of 2-dimensional arrays matching the coordinates of the init array, or a similar 3 dimensional array where the 3rd dimension is the time-step.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.HumanLayer",
    "page": "Dispersal",
    "title": "Dispersal.HumanLayer",
    "category": "type",
    "text": "struct HumanLayer{T} <: Dispersal.AbstractHumanLayer\n\ndata\nAny 2-dimensional AbstractArray matching the coordinates of the init array\n\n\n\n"
},

{
    "location": "index.html#Types-1",
    "page": "Dispersal",
    "title": "Types",
    "category": "section",
    "text": "AbstractLayer \nAbstractSuitabilityLayer \nSuitabilityLayer \nAbstractSuitabilitySequence \nSuitabilitySequence\nHumanLayer"
},

{
    "location": "index.html#Dispersal.suitability",
    "page": "Dispersal",
    "title": "Dispersal.suitability",
    "category": "function",
    "text": "suitability(layers, row::Int, col::Int, t::Number)\n\nReturns a scalar representing cell suitability from one or multiple layers.\n\nFor multiple layers the product of each individual scalar is returned.\n\nLayers of type other than AbstractSuitabilityLayer return 1.0.\n\nArguments\n\nlayers : a single layer or tuple of layers of any type\nrow::Int\ncol::Int\nt::Number : current timestep for interploating layere sequences.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.human_impact",
    "page": "Dispersal",
    "title": "Dispersal.human_impact",
    "category": "function",
    "text": "human_impact(layers, row::Int, col::Int, t)\n\nReturns a scalar indicating human impact from one or multiple layers.\n\nFor multiple layers the product of each individual scalar is returned.\n\nLayers of type other than AbstractHumanLayer return 1.0.\n\nArguments\n\nlayers : a single layer or tuple of layers of any type\nrow::Int\ncol::Int\nt::Number : current timestep for interploating layere sequences.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.sequence_interpolate",
    "page": "Dispersal",
    "title": "Dispersal.sequence_interpolate",
    "category": "function",
    "text": "sequence_interpolate(layer, row, col, t)\n\nInterpolates between layers in a sequence.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.cyclic",
    "page": "Dispersal",
    "title": "Dispersal.cyclic",
    "category": "function",
    "text": "Cycles a time position through a particular timespan length\n\n\n\n"
},

{
    "location": "index.html#Methods-2",
    "page": "Dispersal",
    "title": "Methods",
    "category": "section",
    "text": "suitability\nhuman_impact\nsequence_interpolate\ncyclic"
},

]}
