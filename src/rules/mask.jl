
"""
Simple masking layers. Initialize with a layer of 1s and 0s, 
matching your init matrix dimensions.

`Missing` cannot be used in simulation frames as it will be propagated by 
neighborhood and jump style rule, so active masking is the alternative.
$(FIELDDOCTABLE)
"""
@Layers struct Mask{} <: AbstractCellRule end

@inline applyrule(rule::Mask, data, state, index, args...) = 
    state * get_layers(rule, data, index)
