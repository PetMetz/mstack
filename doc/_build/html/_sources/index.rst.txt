.. MStack documentation master file, created by
   sphinx-quickstart on Sun Apr 02 21:42:34 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*MStack*: stacking disorder tools for *Python*
==================================================
| Peter C. Metz\*
| *Inamori School of Engineering, Alfred University, 2 Pine St., Alfred, NY 14802*
| *\*Contact: pcm1@alfred.edu | (315) 350 1585*


Abstract:
---------------------------------------------------------------
In order to refine stacking disorder models in real and reciprocal space, MSTACK has been written to
extend two established profile generators: DIFFaX, a reciprocal space intensity distribution calculator
built on a stochastic stacking disorder model description; and DiffPy-CMI, a suite of tools including pair
distribution function calculators. MSTACK includes tools to expand the stochastic stacking model parameters
typical of DIFFaX into supercell models suitable for calculation of stacking disordered pair distribution function
data, and to drive refinement of layer structure models from real and reciprocal space data.

MSTACK has been designed with advanced refinement tools in mind. MSTACK is built on the code lmfit which
permits the user to include arbitrary constraint equations enabling parametric refinement. Further, MSTACK
is designed to be compatible with any minimizer method in the SciPy package, enabling the application of
global minimization techniques. Currently implemented minimization methods include the L-BFGS-B non-linear
optimization algorithm and the Differential Evolution algorithm. Finally, lmfit and MSTACK enable the user to
apply Markov-Chain Monte Carlo analysis, via the package emcee, to the resulting fit. This Bayesian statistical
analysis tool has been used predominantly in the astronomical community to interpret data with substantial noise
where the error in the data is uncertain. Application of this tool to PDF data suggests many model parameters are
not normally distributed- an insight that is expected to have substantial impact on the future of scattering analysis
of nanostructured material. 


Acknowldgements:
---------------------------------------------------------------
.. include::
	rst/acknowledgement.rst

Contents:
---------------------------------------------------------------

.. toctree::
	:maxdepth: 3
	:glob: 

	rst/background
	rst/interface
	rst/pairdistributionfunction
	rst/refinement
	rst/structure
	rst/supercell
	rst/transition
	rst/utilities
	

.. background
.. ---------------------------------------------------------------
.. 
.. .. automodule:: mstack.background
.. 	:members:
.. 	:undoc-members:
.. 	:show-inheritance:
.. 	:noindex:
.. 
.. interface
.. ---------------------------------------------------------------
.. 
.. .. automodule:: mstack.interface
.. 	:members:
.. 	:undoc-members:
.. 	:show-inheritance:
.. 	:noindex:
.. 
.. pairdistributionfunction
.. ---------------------------------------------------------------
.. 
.. .. automodule:: mstack.pairdistributionfunction
.. 	:members:
.. 	:undoc-members:
.. 	:show-inheritance:
.. 	:noindex:
.. 
.. refinement
.. ---------------------------------------------------------------
.. 
.. .. automodule:: mstack.refinement
.. 	:members:
.. 	:undoc-members:
.. 	:show-inheritance:
.. 	:noindex:
.. 
.. structure
.. ---------------------------------------------------------------
.. 
.. .. automodule:: mstack.structure
.. 	:members:
.. 	:undoc-members:
.. 	:show-inheritance:
.. 	:noindex:
.. 
.. supercell
.. ---------------------------------------------------------------
.. 
.. .. automodule:: mstack.supercell
.. 	:members:
.. 	:undoc-members:
.. 	:show-inheritance:
.. 	:noindex:
..  
.. transition
.. ---------------------------------------------------------------
.. 
.. .. automodule:: mstack.transition
.. 	:members:
.. 	:undoc-members:
.. 	:show-inheritance:
.. 	:noindex:
.. 
.. utilities
.. ---------------------------------------------------------------
.. 
.. .. automodule:: mstack.utilities
.. 	:members:
.. 	:undoc-members:
.. 	:show-inheritance:
.. 	:noindex:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

