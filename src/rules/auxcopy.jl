"""
    AuxCopy(auxkey, timeindex)

A simple rule that copies aux array slices to a grid over time.
This can be used for comparing simulation dynamics to aux data dynamics.
"""
struct AuxCopy{R,W,L,TI} <: Rule{R,W} 
    "Key for auxillary data source"
    auxkey::L 
    "Precalculated interpolation indices"
    auxtimeindex::TI
end
AuxCopy(args...; kw...) = AuxCopy{Tuple{},:_default_}(args...; kw...)
AuxCopy{R,W}(; auxkey, auxtimeindex=1) where {R,W} = 
    AuxCopy{R,W}(auxkey, auxtimeindex)

DynamicGrids.applyrule(data, rule::AuxCopy, state, cellindex) =
    auxval(data, rule.auxkey, cellindex..., rule.auxtimeindex)

DynamicGrids.precalcrules(rule::AuxCopy, data) =
    precalc_auxtimeindex(aux(data)[rule.auxkey], rule, data)
