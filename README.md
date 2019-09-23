# Dispersal

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cesaraustralia.github.io/Dispersal.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cesaraustralia.github.io/Dispersal.jl/dev)
[![Build Status](https://travis-ci.org/cesaraustralia/Dispersal.jl.svg?branch=master)](https://travis-ci.org/cesaraustralia/Dispersal.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/648h30ifo85wnvfk?svg=true)](https://ci.appveyor.com/project/cesaraustralia/dispersal-jl)
[![Coverage Status](https://coveralls.io/repos/github/cesaraustralia/Dispersal.jl/badge.svg?branch=master)](https://coveralls.io/github/cesaraustralia/Dispersal.jl?branch=master)
[![codecov.io](http://codecov.io/github/cesaraustralia/Dispersal.jl/coverage.svg?branch=master)](http://codecov.io/github/cesaraustralia/Dispersal.jl?branch=master)

[Dispersal.jl](https://github.com/cesaraustralia/Dispersal.jl) extends
[DynamicGrids.jl](https://github.com/cesaraustralia/DynamicGrids.jl) for grid-based
simulations of organism dispersal. Dispersal.jl provides a range of simulation
rules than can be combined to develop complex organism dispersal models. Growth
rates, dispersal kernels, Allee effects, and randomised jump and human assisted
dispersal rules are provided.

These rules can be chained arbitrarily, and custom rules can easily added and
combined with the provided set.

DynamicGrids.jl includes REPLOutput for live display in the REPL, while
[DynamicGridsGtk](https://github.com/cesaraustralia/DynamicGridsGtk.jl) provides
GtkOutput for a simple graphical viewer, and
[DynamicGridsInteract](https://github.com/cesaraustralia/DynamicGridsInteract.jl)
provides an interactive interface for desktop and online web applications.

[GrowthMaps.jl](https://github.com/cesaraustralia/GrowthMaps.jl) can efficiently generate
the layers required for suitability growth rules based on temperature response and stress factors.
