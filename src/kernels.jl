
build_dispersal_kernel(f, params, init, cellsize, r) = begin
    params = typeof(params) <: Tuple ? params : (params,)
    sze = 2r + one(r)
    kernel = zeros(Float64, sze, sze)
    # Paper: l. 97
    for y = -r:r, x = -r:r
        kernel[y+r+one(y), x+r+one(x)] = f(sqrt(y^2 + x^2) * cellsize, params...)
    end
    # Normalise
    kernel ./= sum(kernel)
    typeof(init).name.wrapper(kernel)
end

# Paper: l. 96
exponential(d, a) = exp(-d / a)

