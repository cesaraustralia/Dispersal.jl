# Mixins

@premix @columns struct InstrinsicGrowthRate{GR}
    # Field           | Default  | Flatn | Bounds      | Description
    intrinsicrate::GR | 0.1      | true  | (0.0, 10.0) | "Intrinsic rate of growth"
end

@premix @columns struct CarryCap{CC}
    carrycap::CC      | 100000.0 | true  | (0.0, 1e9)  | "Carrying capacity for each cell. Not currently scaled by area."
end

@premix @columns struct Layers{L,TI}
    layer::L          | nothing  | false | _           | "Data layer"
    timeinterp::TI    | 1        | false | _           | "Precalculated interpolation indices"
end


"""
Extends AbstractCellRule for rules of growth dynamics

For best performance these should be chained with other
AbstractCellRule or following an AbstractNeighborhoodRule.
"""
abstract type AbstractGrowthRule <: AbstractCellRule end


"""
Extends AbstractGrowthRule for growth rules using a heterogenous
growth rate layer. 

[GrowthMaps.jl](http://github.com/cesaraustralia/GrowthMaps.jl)
can produce these growth maps from environmental data.
"""
abstract type AbstractGrowthMapRule <: AbstractGrowthRule end

layer(rule::AbstractGrowthMapRule) = rule.layer 
timeinterp(rule::AbstractGrowthMapRule) = rule.timeinterp

DynamicGrids.precalcrule(rule::AbstractGrowthMapRule, data) = 
    DynamicGrids.precalcrule(layer(rule), rule, data)
DynamicGrids.precalcrule(::AbstractMatrix, rule::AbstractGrowthMapRule, data) = rule
DynamicGrids.precalcrule(::AbstractArray{<:Any,3}, rule::AbstractGrowthMapRule, data) = begin
    @set! rule.timeinterp = precalc_time_interpolation(layer(rule), rule, data)
    rule
end


# Exact growth solutions

"""
Simple fixed exponential growth rate using exact solution.
$(FIELDDOCTABLE)
"""
@InstrinsicGrowthRate struct ExactExponentialGrowth{} <: AbstractGrowthRule end

@inline applyrule(rule::ExactExponentialGrowth, data, state, args...) =
    state * exp(rule.intrinsicrate * timestep(data))

"""
Simple fixed logistic growth rate using exact solution
$(FIELDDOCTABLE)
"""
@CarryCap @InstrinsicGrowthRate struct ExactLogisticGrowth{} <: AbstractGrowthRule end

@inline applyrule(rule::ExactLogisticGrowth, data, state, index, args...) = begin
    @fastmath (state * rule.carrycap) / (state + (rule.carrycap - state) * 
                                         exp(-rule.intrinsicrate * timestep(data)))
end

"""
Exponential growth based on a growth rate layer using exact solution.
$(FIELDDOCTABLE)
"""
@Layers struct ExactExponentialGrowthMap{} <: AbstractGrowthMapRule end

@inline applyrule(rule::ExactExponentialGrowthMap, data, state, index, args...) = begin
    intrinsicrate = layer(rule, data, index)
    @fastmath state * exp(intrinsicrate * timestep(data))
end

"""
Logistic growth based on a growth rate layer, using exact solution.

Saturation only applies with positive growth
$(FIELDDOCTABLE)
"""
@CarryCap @Layers struct ExactLogisticGrowthMap{} <: AbstractGrowthMapRule end

@inline applyrule(rule::ExactLogisticGrowthMap, data, state, index, args...) = begin
    @inbounds intrinsicrate = layer(rule, data, index)
    if intrinsicrate > zero(intrinsicrate)
        @fastmath (state * rule.carrycap) / (state + (rule.carrycap - state) * 
                                             exp(-intrinsicrate * timestep(data)))
    else
        @fastmath state * exp(intrinsicrate * timestep(data))
    end
end

"""
Simple layer mask. Values below a certain threshold are replaced with zero.
$(FIELDDOCTABLE)
"""
@Layers struct MaskGrowthMap{ST} <: AbstractGrowthMapRule
    threshold::ST | 0.5 | true  | (0.0, 1.0) | "Minimum viability index."
end

@inline applyrule(rule::MaskGrowthMap, data, state, index, args...) =
    layer(rule, data, index) >= rule.threshold ? state : zero(state)
