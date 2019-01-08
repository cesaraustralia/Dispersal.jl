
"Simple masking layers. Initialize with a layer of 1s and 0s"
@Layers struct Mask{} <: AbstractCellModel end

@inline rule(model::Mask, data, state, index, args...) = 
    state * get_layers(model.layers, index, data.t)
