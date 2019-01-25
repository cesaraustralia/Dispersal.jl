using SafeTestsets

@time @safetestset "downsampling" begin include("downsampling.jl") end
@time @safetestset "layers" begin include("layers.jl") end
@time @safetestset "human dispersal" begin include("human.jl") end
@time @safetestset "jump dispersal" begin include("jump.jl") end
@time @safetestset "kernels" begin include("kernels.jl") end
@time @safetestset "allee effects" begin include("allee.jl") end
@time @safetestset "growth" begin include("growth.jl") end
@time @safetestset "mask" begin include("mask.jl") end
@time @safetestset "integration" begin include("integration.jl") end
