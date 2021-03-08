
"""
    KernelFormulation

Abstract supertype for functors that calculate the probability density based on distance.

Concrete implementations must define functor methods with the form:

```julia
(k::SomeKernel)(distance) = # do something with `distance` and `k`"
```

Using an anonymous function would not allow rebuildable model parameters.
"""
abstract type KernelFormulation end

"""
    ExponentialKernel <: KernelFormulation

    ExponentialKernel(λ)

Probability density function of distance ``d``.

```math
y = e^{-d/λ}
```

where λ is a shape parameter.
"""
Base.@kwdef struct ExponentialKernel{P} <: KernelFormulation
    λ::P = Param(1.0, bounds=(0.0, 2.0))
end
(f::ExponentialKernel)(d) = exp(-d / f.λ)

"""
    GeometricKernel <: KernelFormulation

    GeometricKernel(α)

Probability density function of distance ``d``.

The Geometric kernel has a power-law decrease.

```math
y = (1+d)^α (α+1)(α+2) / (2 π)
```

where α is a shape parameter.
"""
Base.@kwdef struct GeometricKernel{P} <: KernelFormulation
    α::P = Param(1.0, bounds=(-1000.0, 1000.0))
end
(f::GeometricKernel)(d) = (1 + d)^f.α * ((f.α + 1)*(f.α + 2)) / (2 * π)

"""
    GaussianKernel <: KernelFormulation

    GaussianKernel(α)

Probability density function of distance ``d``.

```math
y = 1/ (π α^2) e^{-d^2/α^2} 
```

where α is a positive parameter.
"""
Base.@kwdef struct GaussianKernel{P} <: KernelFormulation
    α::P = Param(1.0, bounds=(0.0, 1000.0))
end
(f::GaussianKernel)(d) = 1 / (π * f.α^2) * exp(- d^2 / f.α^2)

"""
    WeibullKernel <: KernelFormulation

    WeibullKernel(α,β)

Probability density function of distance ``d``.

```math
y =β /(2 π α^2) d^{β-2} e^{ -d^β/α^β} 
```

where α and β are positive parameters.
"""
Base.@kwdef struct WeibullKernel{A,B} <: KernelFormulation
    α::A = Param(1.0, bounds=(0.0, 1000.0))
    β::B = Param(2.0, bounds=(0.0, 1000.0))
end
(f::WeibullKernel)(d) = f.β / (2 * π * f.α^2) * d^(f.β - 2) * exp(- d^f.β / f.α^f.β)