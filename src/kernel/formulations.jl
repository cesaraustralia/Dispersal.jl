
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
@Base.kwdef struct ExponentialKernel{P} <: KernelFormulation
    "Parameter for adjusting spread of dispersal propability"
    λ::P = Param(1.0, bounds=(0.0, 2.0))
end
(f::ExponentialKernel)(d) = exp(-d / f.λ)

"""
    GeometricKernel <: KernelFormulation

    GeometricKernel(α)

Probability density function of distance ``d``.

The Geometric kernel has a power-law decrease and is isotropic.

```math
y = \frac{(α+1)(α+2)}{2 π} (1+d)^α
```

where α is a shape parameter.
"""
@Base.kwdef struct GeometricKernel{P} <: KernelFormulation
    α::P = Param(1.0, bounds=(-1000.0, 1000.0))
end
(f::GeometricKernel)(d) = (1 + d)^f.α * ((f.α + 1)*(f.α + 2)) / (2 * π)


"""
    GaussianKernel <: KernelFormulation

    GaussianKernel(α)

Probability density function of distance ``d``.

```math
y = \frac{1}{π α^2} \exp\left( \frac{-d^2}{α^2}  \right) 
```

where α is a positive parameter.
"""
@Base.kwdef struct GaussianKernel{P} <: KernelFormulation
    α::P = Param(1.0, bounds=(0.0, 1000.0))
end
(f::GaussianKernel)(d) = 1 / (π * f.α^2) * exp(- d^2 / f.α^2)


"""
    WeibullKernel <: KernelFormulation

    WeibullKernel(α,β)

Probability density function of distance ``d``.

```math
y = \frac{β}{2 π α^2} d^{β-2} \exp\left( \frac{-d^β}{α^β}  \right) 
```

where α and β are positive parameters.
"""
@Base.kwdef struct WeibullKernel{A,B} <: KernelFormulation
    α::A = Param(1.0, bounds=(0.0, 1000.0))
    b::B = Param(1.0, bounds=(0.0, 1000.0))
end
(f::WeibullKernel)(d) = f.β / (2 * π * f.α^2) * d^(f.β - 2) * exp(- d^f.β / f.α^f.b)
