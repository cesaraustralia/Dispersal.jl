
abstract type AbstractDownsampling end


downsample(a::AbstractMatrix, aggregator, scale) = downsample!(initdownsample(a, scale), a, aggregator, scale)

downsample!(down, a::AbstractMatrix, aggregator, scale) = begin
    scale == 1 && return down .= a

    h, w = size(down)
    for i = 1:h, j = 1:w
        i1, j1 = upsample_index((i, j), scale)
        i2, j2 = min.(size(a), (i1, j1) .+ scale .- 1)
        cells = collect(skipmissing(a[i1:i2, j1:j2]))
        down[i, j] = length(cells) > 0 ? aggregator(cells) : missing
    end
    down
end

initdownsample(a, scale) = similar(a, downsample_index(size(a), scale))

upsample_index(index, scale)   = ((index .- 1) .* scale) .+ 1
downsample_index(index, scale) = ((index .- 1) .÷ scale) .+ 1
