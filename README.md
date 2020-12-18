# Dispersal

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cesaraustralia.github.io/Dispersal.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cesaraustralia.github.io/Dispersal.jl/dev)
[![Build Status](https://travis-ci.org/cesaraustralia/Dispersal.jl.svg?branch=master)](https://travis-ci.org/cesaraustralia/Dispersal.jl)
[![codecov.io](http://codecov.io/github/cesaraustralia/Dispersal.jl/coverage.svg?branch=master)](http://codecov.io/github/cesaraustralia/Dispersal.jl?branch=master)

![Spotted wing drosophola dispersal](https://raw.githubusercontent.com/cesaraustralia/packagegifs/master/SWD_USA.gif)

*Fitted dispersal simulation of the spotted-wing drosophola invasion of continental USA.*

[Dispersal.jl](https://github.com/cesaraustralia/Dispersal.jl) extends
[DynamicGrids.jl](https://github.com/cesaraustralia/DynamicGrids.jl) to provide
model components and tools for grid-based simulations of organism dispersal.

Growth rates, dispersal kernels, Allee effects, and randomised jump and human
assisted dispersal rules are provided. These components can be combined into
complex dispersal models. Custom rules can easily added and combined with the
provided set. See the documentation for examples and the lists of included
rules.

There are also methods for optimisation of stochastic dispersal models using
observation data.

[DynamicGridsInteract](https://github.com/cesaraustralia/DynamicGridsInteract.jl)
provides an interactive interface for atom and jupyter notbooks (InteractOuput),
desktop (ElectronOutput) and online web applications (ServerOuput), 
where complete models, including your custom rules, can be manipulated during live
simulations.

[DynamicGridsGtk](https://github.com/cesaraustralia/DynamicGridsGtk.jl) provides
GtkOutput for a simple graphical viewer.

[GrowthMaps.jl](https://github.com/cesaraustralia/GrowthMaps.jl) can efficiently
generate the layers required for growth rules based on temperature
response and stress factors.

