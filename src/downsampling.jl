
"""
    downsample(a::AbstractMatrix, aggregator, scale)

Allocate a new output array and run [`downsample!`](@ref).
"""
function downsample(a::AbstractMatrix, aggregator, scale) 
    downsample!(initdownsample(a, scale), a, aggregator, scale)
end

"""
    downsample!(out::AbstractMatrix, a::AbstractMatrix, aggregator, scale)

Downsample matrix `a` to another matrix `out` of the correct size.

- `aggregator` is a function such as mean or sum that can combine the 
    value of multiple cells to generate the downsampled cell.
- `scale` is the downsampling factor.
"""
function downsample!(out::AbstractMatrix, a::AbstractMatrix, aggregator, scale)
    scale == 1 && return out .= a

    h, w = size(out)
    for i = 1:h, j = 1:w
        i1, j1 = upsample_index((i, j), scale)
        i2, j2 = min.(size(a), (i1, j1) .+ scale .- 1)
        cells = collect(skipmissing(a[i1:i2, j1:j2]))
        out[i, j] = length(cells) > 0 ? aggregator(cells) : missing
    end
    return out
end

"""
    initdownsample(a, scale)

Generate an array for downsampleing array `a` by `scale`.
"""
initdownsample(a, scale) = similar(a, downsample_index(size(a), scale))

"""
    upsample_index(index, scale)

Convert indicies from the downsampled array to the larger original array.
"""
upsample_index(index, scale) = ((index .- 1) .* scale) .+ 1

"""
    downsample_index(index, scale)
Convert indicies from the original array to the downsampled array.
"""
downsample_index(index, scale) = ((index .- 1) .÷ scale) .+ 1
