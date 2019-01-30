using SafeTestsets

@time @safetestset "downsampling" begin include("downsampling.jl") end
@time @safetestset "layers" begin include("layers.jl") end
@time @safetestset "human dispersal" begin include("models/human.jl") end
@time @safetestset "jump dispersal" begin include("models/jump.jl") end
@time @safetestset "kernels" begin include("models/kernels.jl") end
@time @safetestset "allee effects" begin include("models/allee.jl") end
@time @safetestset "growth" begin include("models/growth.jl") end
@time @safetestset "mask" begin include("models/mask.jl") end
@time @safetestset "integration" begin include("integration.jl") end
@time @safetestset "cell optimisation" begin include("optimisation/cell.jl") end
@time @safetestset "region optimisation" begin include("optimisation/region.jl") end
