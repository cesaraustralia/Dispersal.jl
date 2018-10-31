
@columns struct HudginsDispersal{N} <: AbstractOutwardsDispersal
    neighborhood::N = nothing            | true | _
    pop_threshold::Float32 = 0.0006227f0 | true | (0.0, 0.1) 
    growthrate::Float32 = 2.4321f0       | true | (0.0, 10.0)
end

hudgins_precalc(init, suit, human) = begin
    # distances = similar(init, Float32, (size(init) .* 2)...)
    precalc = similar(init, Float32)
    d = similar(init, Float32)
    f = similar(init, Float32)
    zi = similar(init, Float32)
    h, w = convert.(Int32, size(init))
    indices = broadcastable_indices(init)

    broadcast!(zi, suit, human) do s, h
        -0.8438f0 .* s .- 0.1378f0 .* h
    end

    broadcast!(f, zi) do z 
        2.0f0 * 1.1248f0 * CUDAnative.exp(z)/(oneunit(z) + CUDAnative.exp(z))
    end

    for i in one(h):h 
        # println("Cell: ", i)
        for j in one(w):w
            broadcast!(d, indices, (i,), (j,)) do ii, jj, i, j
                CUDAnative.sqrt(convert(Float32, (i - ii) * (i - ii) + (j - jj) * (j - jj)))
            end
            precalc[i, j] = sum(CUDAnative.exp.(-d .* f))
        end
    end
    typeof(zi)(precalc)
end

"""
Dispersal function taken from Hudgins, 'Predicting the spread of all
invasive forest pests in the United States', 2017
"""
neighbors(hood, model::HudginsDispersal, state, indices, t, 
          source, dest, layers, precalc, args...) = begin
    # Ignore cells below the population threshold
    # state > model.pop_threshold || return zero(state)

    # Setup
    # h, w = convert.(Int32, size(precalc))
    # propagules = zero(state)

    # local d::Float32
    # Disperse to the entire grid
    # for i = one(h):h, j = one(w):w
        # Invade the cell
        # @inbounds ZI = -0.8438f0 * suitability(layers, (i, j), t) + -0.1378f0 * human_impact(layers, (i, j), t)
        # f = 2.0f0 * 1.1248f0 * CUDAnative.exp(ZI)/(1.0f0+CUDAnative.exp(ZI))
        # d = CUDAnative.sqrt(convert(Float32, ((row - i)^2 + (col - j)^2)))
        # t = exp(-d * f)/precalc[i, j]

        # invaders = state
        # dest[i, j] += invaders   
        # propagules = propagules + invaders
    # end
    # typeof(state)(propagules)
    state
end
