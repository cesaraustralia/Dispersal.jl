
spec_rand(source, typ, args...) = rand(typ)


abstract type AbstractDownsampling end

struct MeanDownsampling <: AbstractDownsampling end
struct SumDownsampling <: AbstractDownsampling end

aggregate(::MeanDownsampling, xs) = mean(xs)
aggregate(::SumDownsampling, xs) = sum(xs)

downsample(a::AbstractMatrix, mode, scale) = begin
    scale == 1 && return a
    down = similar(a, ((size(a) .- 1) .÷ scale .+ 1))
    h, w = size(down)
    # Add padding columns of missing
    ax = eltype(a)[missing for i = 1:scale * h, j = 1:scale * w] 
    ax[1:size(a, 1), 1:size(a, 2)] .= a
    for i = 1:h, j = 1:w
        i1 = scale * i - 1
        j1 = scale * j - 1
        cells = ax[i1, j1], ax[i1+1, j1], ax[i1, j1+1], ax[i1+1, j1+1]
        if length(collect(skipmissing(cells))) > 0
            down[i, j] = aggregate(mode, skipmissing(cells))
        end
    end
    down
end

upsample_index(index, scale) = ((index .- 1) .* scale) .+ 1
downsample_index(index, scale) = ((index .- 1) .÷ scale) .+ 1
