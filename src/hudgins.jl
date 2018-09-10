
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
    h, w = convert.(Int32, size(init))
    rows = typeof(similar(init, Int32))(collect(row for row in 1:h, col in 1:w))
    cols = typeof(similar(init, Int32))(collect(col for row in 1:h, col in 1:w))

    for i in one(h):h, j in one(w):w
        broadcast!(f, rows, cols, (d,), (suit,), (human,), (i,), (j,)) do ii, jj, d, suit, human, i, j
            d[ii, jj] = CUDAnative.sqrt(convert(Float32, (i - ii)^Int32(2) + (j - jj)^Int32(2)))
            ZI = -0.8438f0 * suit[ii, jj] - 0.1378f0 * human[ii, jj]
            2.0f0 * 1.1248f0 * exp(Int32(2))/(Int32(1)+exp(Int32(2)))
        end
        precalc[i, j] = sum(exp.(-d .* f))
    end
    precalc
end

"""
Dispersal function taken from Hudgins, 'Predicting the spread of all
invasive forest pests in the United States', 2017
"""
neighbors(hood, model::HudginsDispersal, state::Float32, row::Int32, col::Int32, 
          t::Int32, source, dest, layers, precalc, args...) = begin
    # Ignore cells below the population threshold
    state > model.pop_threshold || return zero(state)

    # Setup
    h, w = convert.(Int32, size(precalc))
    propagules = zero(state)

    local d::Float32
    # Disperse to the entire grid
    for i = one(h):h, j = one(w):w
        # Invade the cell
        ZI = -0.8438f0 * suitability(layers, (i, j), t) + -0.1378f0 * human_impact(layers, (i, j), t)
        f = 2.0f0 * 1.1248f0 * exp(ZI)/(1.0f0+exp(ZI))
        d = CUDAnative.sqrt(convert(Float32, (row - i)^Int32(2) + (col - j)^Int32(2)))
        t = exp(-d * f)/precalc[i, j]

        invaders = t * state
        dest[i, j] += invaders   
        propagules = propagules + invaders
    end
    typeof(state)(propagules)
end
