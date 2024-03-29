# Dispersal.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cesaraustralia.github.io/Dispersal.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cesaraustralia.github.io/Dispersal.jl/dev)
![CI](https://github.com/cesaraustralia/Dispersal.jl/workflows/CI/badge.svg)
[![codecov.io](http://codecov.io/github/cesaraustralia/Dispersal.jl/coverage.svg?branch=master)](http://codecov.io/github/cesaraustralia/Dispersal.jl?branch=master)


[Dispersal.jl](https://github.com/cesaraustralia/Dispersal.jl) is a collection
of rules for organism growth and dispersal that can be used for modelling a wide
range of ecological problems, and easily extended and combined with custom rules
and other packages. Dispersal.jl extends
[DynamicGrids.jl](https://github.com/cesaraustralia/DynamicGrids.jl), a
high-performance grid-based modelling framework with clean syntax.



![Spotted wing drosophola dispersal](https://raw.githubusercontent.com/cesaraustralia/packagegifs/master/SWD_USA.gif)

*A simulation of the spotted-wing drosophola invasion of the continental United 
States as detailed in [Maino, Schouten, and Umina (2021)](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/1365-2664.13812)*

Growth rates, dispersal kernels, Allee effects, and randomised jump and human
assisted dispersal rules are provided. These components can be combined into
complex dispersal models. Custom rules can easily added and combined with the
provided set. See the documentation for examples and the lists of included
rules.

[DynamicGridsInteract](https://github.com/cesaraustralia/DynamicGridsInteract.jl)
provides an interactive interface for atom and jupyter notbooks (InteractOuput),
desktop (ElectronOutput) and online web applications (ServerOuput), 
where complete models, including your custom rules, can be manipulated during live
simulations.

[DynamicGridsGtk](https://github.com/cesaraustralia/DynamicGridsGtk.jl) provides
GtkOutput for a simple graphical viewer.

[GrowthMaps.jl](https://github.com/cesaraustralia/GrowthMaps.jl) can efficiently 
generate summarised raster data for vital rates (e.g. intrinsic growth rates) based on 
higher resolution and shorter interval environmental data.  

