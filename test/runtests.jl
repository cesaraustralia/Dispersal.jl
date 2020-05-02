using SafeTestsets

@time @safetestset "optimisation" begin include("optimisation/optimisation.jl") end
@time @safetestset "ouput" begin include("optimisation/output.jl") end
@time @safetestset "integration" begin include("integration.jl") end
@time @safetestset "downsampling" begin include("downsampling.jl") end
@time @safetestset "layers" begin include("layers.jl") end
@time @safetestset "human dispersal" begin include("rules/human.jl") end
@time @safetestset "jump dispersal" begin include("rules/jump.jl") end
@time @safetestset "kernels" begin include("rules/kernels.jl") end
@time @safetestset "allee effects" begin include("rules/allee.jl") end
@time @safetestset "growth" begin include("rules/growth.jl") end
