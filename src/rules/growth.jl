# Mixins

@mix @columns struct CarryCap{CC}
    # Field           | Default  | Flat  | Bounds      | Description
    carrycap::CC      | 100000.0 | true  | (0.0, 1e9)  | "Carrying capacity for each cell. Not currently scaled by area."
end

@mix @Timestep struct InstrinsicGrowthRate{GR}
    intrinsicrate::GR | 0.1      | true  | (0.0, 10.0) | "Intrinsic rate of growth per timestep"
end

"""
Extends CellRule for rules of growth dynamics

For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""
abstract type GrowthRule{R,W} <: CellRule{R,W} end


"""
Extends GrowthRule for growth rules using a heterogenous
growth rate layer.

[GrowthMaps.jl](http://github.com/cesaraustralia/GrowthMaps.jl)
can produce these growth maps from environmental data.
"""
abstract type GrowthMapRule{R,W} <: GrowthRule{R,W} end

# Base.key(rule::GrowthMapRule{R,W,K}) where {R,W,K} = K
# layer(rule::GrowthMapRule, data) = axkey(rule)

DynamicGrids.precalcrules(rule::GrowthMapRule, data) = begin
    if :timestep in fieldnames(typeof(rule))
        rule = precalctimestep(rule, data)
    end
    precalclayer(layer(rule, data), rule, data)
end

DynamicGrids.precalcrules(rule::GrowthRule, data) =
    :timestep in fieldnames(typeof(rule)) ? precalctimestep(rule, data) : rule


# Exact growth solutions

"""
Simple fixed exponential growth rate using exact solution.
$(FIELDDOCTABLE)
"""
@InstrinsicGrowthRate struct ExactExponentialGrowth{R,W} <: GrowthRule{R,W} end

@inline applyrule(data, rule::ExactExponentialGrowth, population, args...) = begin
    population > zero(population) || return zero(population)
    population * exp(rule.intrinsicrate * rule.nsteps)
end

"""
Simple discrete growth rate.
$(FIELDDOCTABLE)
"""

@InstrinsicGrowthRate struct DiscreteGrowth{R,W} <: GrowthRule{R,W} end

@inline applyrule(data, rule::DiscreteGrowth, population, args...) = begin
    population > zero(population) || return zero(population)
    rule.intrinsicrate * population
end

@mix @Timestep struct InstrinsicGrowthRate2{GR1, GR2}
    intrinsicrate1::GR1 | 0.1      | true  | (0.0, 10.0) | "Intrinsic rate of growth per timestep"
    intrinsicrate2::GR2 | 0.1      | true  | (0.0, 10.0) | "Intrinsic rate of growth per timestep"
end

@InstrinsicGrowthRate2  @CarryCap struct DiscreteGrowth2{R,W} <: GrowthRule{R,W} end

@inline applyrule(data, rule::DiscreteGrowth2, (population1, population2), args...) = begin
    ( rule.intrinsicrate1 * population1 * (1-(population1 + population2)/rule.carrycap),
      rule.intrinsicrate2 * population2 * (1-(population1 + population2)/rule.carrycap) )
end

@mix @Timestep struct InstrinsicGrowthRate3{GR1, GR2, GR3}
    intrinsicrate1::GR1 | 0.1      | true  | (0.0, 10.0) | "Intrinsic rate of growth per timestep"
    intrinsicrate2::GR2 | 0.1      | true  | (0.0, 10.0) | "Intrinsic rate of growth per timestep"
    intrinsicrate3::GR3 | 0.1      | true  | (0.0, 10.0) | "Intrinsic rate of growth per timestep"
end


@InstrinsicGrowthRate3  @CarryCap struct DiscreteGrowth3{R,W} <: GrowthRule{R,W} end

@inline applyrule(data, rule::DiscreteGrowth3, (population1,population2,population3), args...) = begin
    carrycapEffect = (1-sumPop/rule.carrycap)
    ( rule.intrinsicrate1 * population1 * carrycapEffect,
      rule.intrinsicrate2 * population2 * carrycapEffect,
      rule.intrinsicrate3 * population3 * carrycapEffect )
end


"""
Simple fixed logistic growth rate using exact solution
$(FIELDDOCTABLE)
"""
@InstrinsicGrowthRate @CarryCap struct ExactLogisticGrowth{R,W} <: GrowthRule{R,W} end

@inline applyrule(data, rule::ExactLogisticGrowth, population, index, args...) = begin
    population > zero(population) || return zero(population)
    @fastmath (population * rule.carrycap) / (population + (rule.carrycap - population) *
                                         exp(-rule.intrinsicrate * rule.nsteps))
end


"""
Simple discrete growth rate.
$(FIELDDOCTABLE)
"""
@Layers struct DiscreteGrowthMap{R,W} <: GrowthMapRule{R,W} end

@inline applyrule(data, rule::DiscreteGrowthMap, population, index, args...) = begin
    population > zero(population) || return zero(population)
    intrinsicrate = layer(rule, data, index)
    @fastmath intrinsicrate * population
end



"""
Exponential growth based on a growth rate layer using exact solution.
$(FIELDDOCTABLE)
"""
@Timestep @Layers struct ExactExponentialGrowthMap{R,W} <: GrowthMapRule{R,W} end

@inline applyrule(data, rule::ExactExponentialGrowthMap, population, index, args...) = begin
    population > zero(population) || return zero(population)
    intrinsicrate = layer(rule, data, index)
    @fastmath population * exp(intrinsicrate * rule.nsteps)
end

"""
Logistic growth based on a growth rate layer, using exact solution.

Saturation only applies with positive growth
$(FIELDDOCTABLE)
"""
@Layers @Timestep @CarryCap struct ExactLogisticGrowthMap{R,W} <: GrowthMapRule{R,W} end

@inline applyrule(data, rule::ExactLogisticGrowthMap, population, index, args...) = begin
    population > zero(population) || return zero(population)
    @inbounds intrinsicrate = layer(rule, data, index)
    if intrinsicrate > zero(intrinsicrate)
        @fastmath (population * rule.carrycap) / (population + (rule.carrycap - population) *
                                             exp(-intrinsicrate * rule.nsteps))
    else
        @fastmath population * exp(intrinsicrate * rule.nsteps)
    end
end

"""
Simple layer mask. Values below a certain threshold are replaced with zero.
$(FIELDDOCTABLE)
"""
@Layers @Timestep struct MaskGrowthMap{R,W,ST} <: GrowthMapRule{R,W}
    threshold::ST | 0.5 | true  | (0.0, 1.0) | "Minimum viability index."
end

@inline applyrule(data, rule::MaskGrowthMap, population, index, args...) =
    layer(rule, data, index) >= rule.threshold ? population : zero(population)


"""
Logistic growth based on a growth rate layer, using exact solution.

$(FIELDDOCTABLE)
"""

@Layers @Timestep struct ExactLogisticGrowthMap2{R,W} <: GrowthMapRule{R,W} end

@inline applyrule(data, rule::ExactLogisticGrowthMap2, population, timeindex, args...) = begin
    population > zero(population) || return zero(population)
    intrinsicrate = layer(rule, data, timeindex, 1)
    carrycap = layer(rule, data, timeindex, 2)
    if intrinsicrate > zero(intrinsicrate)
        @fastmath (population * carrycap) / (population + (carrycap - population) *
                                             exp(-intrinsicrate * rule.nsteps))
    else
        @fastmath population * exp(intrinsicrate * rule.nsteps)
    end
end


# Define the struct from scratch without all the macros, just Base@kwdef for default values
# instead of `layerkey` we have two keys `ratekey` and `carrycapkey` so there are two
# datasets you can retrieve. 
@mix @columns struct DoubleLayers{RK,CK,TI}
    ratekey::RK       | nothing | false | _ | "key for aux layer 1"
    carrycapkey::CK   | nothing | false | _ | "key for aux layer 2"
    timeindex::TI     | 1       | false | _ | "Precalculated interpolation indices. Not set by users"
end

@DoubleLayers @Timestep struct ExactLogisticGrowthMap3{R,W} <: GrowthMapRule{R,W} end
#    ratekey::RK         = Val{:intrinsicrate}
#    carrycapkey::CK     = Val{:carrycap}
#    timeindex::TI       = 1
#    nsteps::N           = 1
# end

@inline applyrule(data, rule::ExactLogisticGrowthMap3, population, index, args...) = begin
    population > zero(population) || return zero(population)
    # I'll think of a cleaner way to do this part
    intrinsicrate = aux(data)[unwrap(rule.ratekey)][index..., timeindex(rule)]
    carrycap = aux(data)[unwrap(rule.carrycapkey)][index..., timeindex(rule)]

    if intrinsicrate > zero(intrinsicrate)
        @fastmath (population * carrycap) / (population + (carrycap - population) *
                                             exp(-intrinsicrate * rule.nsteps))
    else
        @fastmath population * exp(intrinsicrate * rule.nsteps)
    end
end