
"""
Simple masking layers. Initialize with a layer of 1s and 0s, 
matching your init matrix dimensions.

`Missing` cannot be used in simulation frames as it will be propagated by 
neighborhood and jump style models, so active masking is the alternative.
$(FIELDDOCTABLE)
"""
@Layers struct Mask{} <: AbstractCellModel end

@inline rule(model::Mask, data, state, index, args...) = 
    state * get_layers(model, data, index)
