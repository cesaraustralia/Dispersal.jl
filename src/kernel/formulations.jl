
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
