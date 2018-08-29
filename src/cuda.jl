using CuArrays,
      GPUArrays,
      CUDAnative

pressure(model, source, cc, randomstate, args...) = begin
    rnd = GPUArrays.gpu_rand(Float64, CuArrays.CuKernelState(), randomstate)
    CUDAnative.pow(rnd, model.prob_threshold) > (1 - cc) / 1
end
