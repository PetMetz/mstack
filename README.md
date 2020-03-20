# mstack
Not fit for public consumption. Currently refactoring. Please disregard.

I suspect the best option for these tools is for the DIFFaX subprocess to be rewritten as a CMI-compatible Calculator,
which could leverage a number of existing optimizations in the diffpy srfit package. I.e. the pattern doesn't need
to be reevaluated if only a scale factor is modified.

Likewise, the real space stacking tools should become a fit contribution constructor that manages pointers / proxies
for hierarchically constrained variables.

The current version is recently futurized for py3.7 compatibility, which is quite useful as subprocess now has builtin
timeout for spawned processes, and DIFFaX will hang if fed inappropriate input.

P.S. This was the project in which I learned to code in Python. If my code looks like a catastrophe, that's because it is ;)

# mstack

Stacking disorder tools for Python extending the [DiffPy-CMI](http://www.DiffPy.org/products/diffpycmi) 
and [DIFFaX](http://www.public.asu.edu/~mtreacy/DIFFaX.html) profile generators.

# Abstract

In order to refine stacking disorder models in real and reciprocal space, MSTACK 
has been written to extend two established profile generators: DIFFaX, a reciprocal
space intensity distribution calculator built on a stochastic stacking disorder model
description; and DiffPy-CMI, a suite of tools including pair distribution function 
calculators. 

MSTACK includes tools to expand the stochastic stacking model parameters typical
of DIFFaX into supercell models suitable for calculation of stacking disordered pair 
distribution function data, and to drive refinement of layer structure models from real
and reciprocal space data. 

MSTACK has been designed with advanced refinement tools in mind. MSTACK is built on
the code lmfit which permits the user to include arbitrary constraint equations enabling
parametric refinement. Further, MSTACK is designed to be compatible with any minimizer
method in the SciPy package, enabling the application of global minimization techniques. 

Currently implemented minimization methods include the L-BFGS-B non-linear optimization
algorithm and the Differential Evolution algorithm. Finally, lmfit and MSTACK enable the user
to apply Markov-Chain Monte Carlo analysis, via the package emcee, to the resulting fit.

# etc.

API has been compiled in .PDF and HTML formats.

Please see the accompanying manual for further information.