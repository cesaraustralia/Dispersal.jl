
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
(f::ExponentialKernel)(distance) = exp(-distance / f.λ)

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
(f::GeometricKernel)(distance) = (1+distance)^f.α * ((f.α+1)*(f.α+2))/(2*π)


"""
    GaussianKernel <: KernelFormulation

    GaussianKernel(α)

Probability density function of distance ``d``.

```math
y = \frac{(α+1)(α+2)}{2 π} (1+d)^α
```

where α is a positive parameter.
"""
@Base.kwdef struct GaussianKernel{P} <: KernelFormulation
    α::P = Param(1.0, bounds=(0.0, 1000.0))
end
(f::GaussianKernel)(d) = 1/(π*f.α^2)*exp(-d^2/f.α^2)


"""
    BivariateStudentKernel <: KernelFormulation

    BivariateStudentKernel(α, β, θ, c1, c2)

Probability density function of distance ``d``.
The Bivariate Student kernel has a power-law decrease and is anisotropic.

# Keyword Arguments
- `α` and `β` are positive parameters.
- `θ` is the angle made by the vector.
- `c1` and `c2` are positive parameter with ``c2 ∈ [0, 2π)``

```math
y = (β-1) / (π α^2) (1+ d^2/α^2) )^(-β) e^(c1 cos(θ-c_2))
```
"""
@Base.kwdef struct BivariateStudentKernel{P1,P2,P3,P4,P5} <: KernelFormulation
    α::P1 = Param(1.0, bounds=(0.0, 1000.0))
    β::P2 = Param(1.0, bounds=(0.0, 1000.0))
    θ::P3 = Param(0.5, bounds=(0.0, 10.0))
    c1::P4 = Param(0.5, bounds=(0.0, 10.0))
    c2::P5 = Param(0.5, bounds=(0.0, 2*π))
end
(f::BivariateStudentKernel)(d) = (f.β-1)/(π*f.α^2)*(1+d^2/f.α^2)*exp(f.c1 * cos(f.θ-f.c2)) 