
"""
Simple discrete growth rate.
$(FIELDDOCTABLE)
"""

struct DiscreteGrowth{R,W,GR} <: GrowthRule{R,W} 
    "Intrinsic rate of growth per timestep"
    intrinsicrate::GR
end
function DiscreteGrowth{R,W}(; intrinsicrate=Param(0.1, bounds=(0.0, 10.0))) where {R,W}
    DiscreteGrowth{R,W}(intrinsicrate)
end

@inline function applyrule(data, rule::DiscreteGrowth, population, args...)
    population > zero(population) || return zero(population)
    rule.intrinsicrate * population
end

struct CappedDiscreteGrowth{R,W,GR} <: GrowthRule{R,W} 
    "Intrinsic rate of growth per timestep"
    intrinsicrate::GR
end
function DiscreteGrowth{R,W}(; intrinsicrate=Param(0.1, bounds=(0.0, 10.0))) where {R,W}
    DiscreteGrowth{R,W}(intrinsicrate)
end

@inline function applyrule(data, rule::CappedDiscreteGrowth, pops, args...)
    carrycap_effect = (1 - sum(pops)/ rule.carrycap)
    map(pops,  rule.intrinsicrate1) do pop, ir
        ir * pop * carrycap_effect
    end
end


"""
Simple discrete growth rate.
$(FIELDDOCTABLE)
"""
# @Layers 
struct DiscreteGrowthMap{R,W} <: GrowthMapRule{R,W} 
end

@inline function applyrule(data, rule::DiscreteGrowthMap, population, index, args...)
    population > zero(population) || return zero(population)
    intrinsicrate = layer(rule, data, index)
    @fastmath intrinsicrate * population
end


# Define the struct from scratch without all the macros, just Base@kwdef for default values
# instead of `layerkey` we have two keys `ratekey` and `carrycapkey` so there are two
# datasets you can retrieve. 
struct DoubleLayers{RK,CK,TI}
end

# @DoubleLayers @Timestep 
struct ExactLogisticGrowthMap3{R,W} <: GrowthMapRule{R,W}
    "Key for growth rate layer"
    ratekey::RK  
    "key for carrycap layer"
    carrycapkey::CK
    "Precalculated interpolation indices. Not set by users"
    timeindex::TI
end
function ExactLogisticGrowthMap3{R,W}(
    ratekey=Val{:intrinsicrate},
    carrycapkey=Val{:carrycap},
    timeindex=1,
    nsteps=1
) where {R,W}
    ExactLogisticGrowthMap3{R,W}(ratekey, carrycapkey, timeindex, nsteps)
end

@inline function applyrule(data, rule::ExactLogisticGrowthMap3, population, index, args...)
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
