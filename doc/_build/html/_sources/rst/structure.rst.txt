.. structure

structure
---------------------------------------------------------------
A structure contains atoms (x ,y, z, occ, type, ADP), and 
a structure contains a unit cell (a, b, c, α, β, γ) and asymmetric
unit (list of atoms).

A phase contains layer structures as well as transitions to
propagate the layer motif in along the perpendicular vector.
The stacking vector is always taken as parallel to **c**.

.. automodule::  structure
	:members: build_cif
	:show-inheritance:

structure.Atom
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass::  structure.Atom
	:members:
	:undoc-members:
	:show-inheritance:

structure.Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass::  structure.Structure
	:members:
	:undoc-members:
	:show-inheritance:

structure.Phase
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass::  structure.Phase
	:members:
	:undoc-members:
	:show-inheritance:



