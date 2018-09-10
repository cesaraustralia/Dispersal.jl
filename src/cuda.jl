using CuArrays,
      GPUArrays,
      CUDAnative

pressure(model, source, cc, randomstate, args...) = begin
    spec_rand(source, Float64, randomstate)
    CUDAnative.pow(rnd, model.prob_threshold) > (one(cc) - cc) / one(cc)
end


spec_rand(source, typ, randomstate, args...) = GPUArrays.gpu_rand(typ, CuArrays.CuKernelState(), randomstate)
