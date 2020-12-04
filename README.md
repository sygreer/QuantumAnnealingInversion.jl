QuantumAnnealingInversion: Using a quantum annlear to solve inverse problems in Julia
===============================

Description
-----------

QuantumAnnealingInversion is a [Julia](http://julialang.org/) module for solving inverse problems. It uses the [ThreeQ.jl package ](https://github.com/omalled/ThreeQ.jl/) to interface with [D-Wave](https://www.dwavesys.com/) hardware. 

A number of [examples](https://github.com/sygreer/QuantumAnnealingInversion.jl/tree/master/examples) are including that illustrate how to use QuantumAnnealingInversion. The two main examples are [solving a 2D hydrology inverse problem](https://github.com/sygreer/QuantumAnnealingInversion.jl/blob/master/examples/hydrology2d/hydrolin5.jl), and [solving a 2D seismic inverse problem](https://github.com/sygreer/QuantumAnnealingInversion.jl/blob/master/examples/seismic2d/seismic.jl).

Installation
------------

First, Julia  must be installed. Julia binaries can be obtained [here](http://julialang.org/downloads/). In order to actually use the D-Wave, D-Wave hardware must be available and at least part of D-Wave's software stack must also be installed. D-Wave's [qbsolv](https://github.com/dwavesystems/qbsolv) tool can be used without D-Wave hardware. Other components (the Python SAPI and dw) must be obtained from D-Wave, if needed.

Second, this repo must be cloned into a place where Julia can find it. You can tell Julia where to look for it by adding a line to your .juliarc.jl file like
```
push!(LOAD_PATH, "/dir/where/this/repo/is")
```
For instance, my .juliarc.jl file contains
```
push!(LOAD_PATH, "$(homedir())/codes")
```
and within ~/codes, there is a directory containing this repo called "QuantumAnnealingInversion.jl" or "QuantumAnnealingInversion".

License
-------

QuantumAnnealingInversion is provided under a BSD-ish license with a "modifications must be indicated" clause.  See LICENSE.md file for the full text.

Authors
------

Sarah Greer, <sygreer@mit.edu>
Daniel O'Malley, <omalled@lanl.gov>

