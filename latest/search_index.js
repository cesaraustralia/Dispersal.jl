var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Dispersal",
    "page": "Introduction",
    "title": "Dispersal",
    "category": "module",
    "text": "This package extends Cellular.jl\n\nDispersal.jl provides a range of dispersal modules that can be combined to build grid-based organism dispersal simulations.\n\nThe design provides a solid framewhork while allowing customisation of any  aspect of a simulation. A model may start with the defaults and formulations provided,  but incrementally customise them for a particular use-case, to any level of detail. \n\nAdditionally, modules, outputs, neighborhoods provided by Cellular.jl or other  packages that extend it may all be incorporated into a simulaiton.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.jl-1",
    "page": "Introduction",
    "title": "Dispersal.jl",
    "category": "section",
    "text": "CurrentModule = DispersalDispersal"
},

{
    "location": "index.html#Example-1",
    "page": "Introduction",
    "title": "Example",
    "category": "section",
    "text": "This is the general pattern of running model in Dispersal.jl. We declare a suitability layer and initial conditions. Then we run a dispersal simulation that combines local and jump dispersal, for three timesteps.using Cellular\nusing Dispersal\n\n# add a 0.0-1.0 scaled raster array that represents habitat suitability\nsuitability = [0.5 0.0 0.3 0.0 0.0 0.3 0.0;\n               0.0 0.2 0.8 0.0 0.9 0.6 0.0;\n               0.0 0.5 1.0 1.0 1.0 0.0 0.0;\n               0.0 0.0 1.0 1.0 0.7 0.0 0.0;\n               0.0 0.0 0.6 1.0 0.0 0.0 0.0;\n               0.0 0.1 0.6 0.0 0.0 0.0 0.0;\n               0.0 0.0 0.0 0.0 0.0 0.0 0.7]\n\n# Make an init array the same size as your suitability layer, and seed it\ninit = [0  0  0  0  0  0  0;\n        0  0  0  0  0  0  0;\n        0  0  0  0  0  0  0;\n        0  0  0  0  0  0  0;\n        0  0  1  0  0  0  0;\n        0  0  0  0  0  0  0;\n        0  0  0  0  0  0  0]\n\n# You also could do this with:\n# init = zeros(suitability); init[5, 3] = 1\n\n\n# Define a dispersal kernel function\n\nf = d -> e^-d\n\n# Define the neighborhood, using the dispersal kernel and a radius\nhood = DispersalNeighborhood(; f=f, radius=2)\n\n# Define additional raster layers\nlayers = SuitabilityLayer(suitability)\n\n# Define disersal modules\nlocaldisp = LocalDispersal(layers=layers, neighborhood=hood)\njumpdisp = JumpDispersal(layers=layers)\n\n# Set the output type\noutput = ArrayOutput(init)\n\n# Run the simulation\nsim!(output, (localdisp, jumpdisp), init; time=1:3) \n\noutput.frames[3]"
},

{
    "location": "index.html#Dispersal.AbstractLocalDispersal",
    "page": "Introduction",
    "title": "Dispersal.AbstractLocalDispersal",
    "category": "type",
    "text": "abstract AbstractLocalDispersal <: Cellular.AbstractModel\n\nExtend to modify LocalDispersal\n\n\n\n"
},

{
    "location": "index.html#Dispersal.LocalDispersal",
    "page": "Introduction",
    "title": "Dispersal.LocalDispersal",
    "category": "type",
    "text": "struct LocalDispersal{N, L, P, S} <: Dispersal.AbstractLocalDispersal\n\nLocal dispersal within a DispersalNeighborhood] or other AbstractNeighborhood\n\nKeyword Arguments\n\nneighborhood: = DispersalNeighborhood()\nprob: A real number between one and zero.\nlayers: AbstractLayers or a single AbstractLayer. The default is nothing.\ncellsize::S = A number or Unitful.jl distance.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.AbstractJumpDispersal",
    "page": "Introduction",
    "title": "Dispersal.AbstractJumpDispersal",
    "category": "type",
    "text": "abstract AbstractJumpDispersal <: Cellular.AbstractInPlaceModel\n\nExtend to modify JumpDispersal\n\n\n\n"
},

{
    "location": "index.html#Dispersal.JumpDispersal",
    "page": "Introduction",
    "title": "Dispersal.JumpDispersal",
    "category": "type",
    "text": "struct JumpDispersal{L, S} <: Dispersal.AbstractJumpDispersal\n\nLocal dispersal within a DispersalNeighborhood] or other AbstractNeighborhood\n\nKeyword Arguments\n\nspotrange: A number or Unitful.jl distance.\nprob: A real number between one and zero.\nlayers: AbstractLayers or a single AbstractLayer. The default is nothing.\ncellsize::S = A number or Unitful.jl distance.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.AbstractHumanDispersal",
    "page": "Introduction",
    "title": "Dispersal.AbstractHumanDispersal",
    "category": "type",
    "text": "abstract AbstractHumanDispersal <: Cellular.AbstractInPlaceModel\n\nInherit to exten human dispersal\n\n\n\n"
},

{
    "location": "index.html#Dispersal.HumanDispersal",
    "page": "Introduction",
    "title": "Dispersal.HumanDispersal",
    "category": "type",
    "text": "struct HumanDispersal{L, S} <: Dispersal.AbstractHumanDispersal\n\nHuman dispersal model.\n\nhuman\nlayers\ncellsize\n\n\n\n"
},

{
    "location": "index.html#Cellular.rule",
    "page": "Introduction",
    "title": "Cellular.rule",
    "category": "function",
    "text": "rule(model, state, index, t, args)\n\n\nShort range rule for AbstractLocalDispersal dispersal. Cells are invaded  if there is pressure and suitable habitat, otherwise left as-is.\n\n\n\nrule(model, state, index, t, source, args)\n\n\nLong range rule for AbstractJumpDispersal. A random cell within the spotrange is invaded if it is suitable.\n\n\n\n"
},

{
    "location": "index.html#Models-1",
    "page": "Introduction",
    "title": "Models",
    "category": "section",
    "text": "AbstractLocalDispersal\nLocalDispersal\nAbstractJumpDispersal\nJumpDispersal\nAbstractHumanDispersal\nHumanDispersalThese modules all trigger methods of rule:rule"
},

{
    "location": "index.html#Dispersal.AbstractDispersalNeighborhood",
    "page": "Introduction",
    "title": "Dispersal.AbstractDispersalNeighborhood",
    "category": "type",
    "text": "abstract AbstractDispersalNeighborhood <: Cellular.AbstractNeighborhood\n\nNeighborhoods for dispersal\n\n\n\n"
},

{
    "location": "index.html#Dispersal.DispersalNeighborhood",
    "page": "Introduction",
    "title": "Dispersal.DispersalNeighborhood",
    "category": "type",
    "text": "struct DispersalNeighborhood{K, S} <: Dispersal.AbstractDispersalNeighborhood\n\nA neighborhood built from a dispersal kernel function and a cell radius.\n\nArguments:\n\ndispkernel::AbstractArray{T,2}\noverflow::AbstractOverflow\n\n\n\n"
},

{
    "location": "index.html#Neighborhoods-1",
    "page": "Introduction",
    "title": "Neighborhoods",
    "category": "section",
    "text": "Extend Cellular.AbstractNeighborhoodAbstractDispersalNeighborhood\nDispersalNeighborhood"
},

{
    "location": "index.html#Cellular.neighbors",
    "page": "Introduction",
    "title": "Cellular.neighbors",
    "category": "function",
    "text": "neighbors(h, state, index, t, source, args)\n\n\nReturns nieghbors for a DispersalNeighborhood, looping over the array of dispersal propabilities.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.build_dispersal_kernel",
    "page": "Introduction",
    "title": "Dispersal.build_dispersal_kernel",
    "category": "function",
    "text": "build_dispersal_kernel(f, r)\n\n\nAccepts a dispersal kernel function and integer radius, and returns an array of probabilities, of size  2r + 1 * 2r + 1.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.pressure",
    "page": "Introduction",
    "title": "Dispersal.pressure",
    "category": "function",
    "text": "pressure(model, cc)\n\ndefined at Dispersal/src/disperse.jl:122.\n\nCalculates the propagule pressure from the output of a neighborhood.\n\n\n\n"
},

{
    "location": "index.html#Customisation-1",
    "page": "Introduction",
    "title": "Customisation",
    "category": "section",
    "text": "neighbors\nbuild_dispersal_kernel\npressure"
},

{
    "location": "index.html#Dispersal.AbstractLayer",
    "page": "Introduction",
    "title": "Dispersal.AbstractLayer",
    "category": "type",
    "text": "abstract AbstractLayer\n\nLayers wrap an array, normally built from some kind of raster.\n\nThe wrapper defines the purpose of the layer and allows specialised  method dispatch to utilise them.\n\n\n\n"
},

{
    "location": "index.html#Dispersal.AbstractSuitabilityLayer",
    "page": "Introduction",
    "title": "Dispersal.AbstractSuitabilityLayer",
    "category": "type",
    "text": "abstract AbstractSuitabilityLayer <: Dispersal.AbstractLayer\n\n\n\n"
},

{
    "location": "index.html#Dispersal.SuitabilityLayer",
    "page": "Introduction",
    "title": "Dispersal.SuitabilityLayer",
    "category": "type",
    "text": "struct SuitabilityLayer{T} <: Dispersal.AbstractSuitabilityLayer\n\ndata\n\n\n\n"
},

{
    "location": "index.html#Dispersal.AbstractSuitabilitySequence",
    "page": "Introduction",
    "title": "Dispersal.AbstractSuitabilitySequence",
    "category": "type",
    "text": "abstract AbstractSuitabilitySequence <: Dispersal.AbstractSuitabilityLayer\n\n\n\n"
},

{
    "location": "index.html#Dispersal.SuitabilitySequence",
    "page": "Introduction",
    "title": "Dispersal.SuitabilitySequence",
    "category": "type",
    "text": "struct SuitabilitySequence{T, D} <: Dispersal.AbstractSuitabilitySequence\n\ntimespan\ndata\n\n\n\n"
},

{
    "location": "index.html#Dispersal.HumanLayer",
    "page": "Introduction",
    "title": "Dispersal.HumanLayer",
    "category": "type",
    "text": "struct HumanLayer{T} <: Dispersal.AbstractLayer\n\ndata\n\n\n\n"
},

{
    "location": "index.html#Dispersal.AbstractLayers",
    "page": "Introduction",
    "title": "Dispersal.AbstractLayers",
    "category": "type",
    "text": "abstract AbstractLayers\n\n\n\n"
},

{
    "location": "index.html#Dispersal.Layers",
    "page": "Introduction",
    "title": "Dispersal.Layers",
    "category": "type",
    "text": "struct Layers{S, F} <: Dispersal.AbstractLayers\n\nsuitability\nhuman\n\n\n\n"
},

{
    "location": "index.html#Layers-1",
    "page": "Introduction",
    "title": "Layers",
    "category": "section",
    "text": "Layers are a concept not present in Cellular.jl. They provide  overlay grids of additional information about dispersal potential. Raw arrays are wrapped in a type that determines how they will be used.AbstractLayer \nAbstractSuitabilityLayer \nSuitabilityLayer \nAbstractSuitabilitySequence \nSuitabilitySequence\nHumanLayer\nAbstractLayers\nLayers"
},

{
    "location": "index.html#Dispersal.suitability",
    "page": "Introduction",
    "title": "Dispersal.suitability",
    "category": "function",
    "text": "suitability(layers, row, col, t)\n\n\nReturns a suitability scalar from a single layer, or the product of multiple layers\n\n\n\n"
},

{
    "location": "index.html#Dispersal.seq_interpolate",
    "page": "Introduction",
    "title": "Dispersal.seq_interpolate",
    "category": "function",
    "text": "seq_interpolate(layer, row, col, t)\n\n\nInterpolates between layers in a sequence\n\n\n\n"
},

{
    "location": "index.html#Dispersal.get_cell",
    "page": "Introduction",
    "title": "Dispersal.get_cell",
    "category": "function",
    "text": "get_cell(layer, row, col, pos)\n\n\nReturn a particular cell from a layer, given row, column and timestep)\n\n\n\n"
},

{
    "location": "index.html#Customisation-2",
    "page": "Introduction",
    "title": "Customisation",
    "category": "section",
    "text": "suitability\nseq_interpolate\nget_cell"
},

]}
