using Dispersal, Aqua, SafeTestsets

if VERSION >= v"1.5.0"
    # Amibiguities are not owned by Dispersal
    # Aqua.test_ambiguities([Dispersal, Base, Core])
    Aqua.test_unbound_args(Dispersal)
    Aqua.test_undefined_exports(Dispersal)
    Aqua.test_project_extras(Dispersal)
    Aqua.test_stale_deps(Dispersal)
    Aqua.test_deps_compat(Dispersal)
    Aqua.test_project_toml_formatting(Dispersal)
end

@time @safetestset "integration" begin include("integration.jl") end
@time @safetestset "downsampling" begin include("downsampling.jl") end
@time @safetestset "human dispersal" begin include("human.jl") end
@time @safetestset "jump dispersal" begin include("jump.jl") end
@time @safetestset "kernels" begin include("kernels.jl") end
@time @safetestset "allee effects" begin include("allee.jl") end
@time @safetestset "growth" begin include("growth.jl") end
@time @safetestset "mortality" begin include("mortality.jl") end
@time @safetestset "chain" begin include("chain.jl") end
