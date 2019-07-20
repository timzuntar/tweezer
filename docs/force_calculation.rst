Force calculation
=================

The ``force_calc`` module is used to actually compute forces the traps exert on nearby particles. Assuming all trap parameters have been determined, the procedure is straightforward; trap-particle displacements merely need to be multiplied with trap coefficients, and the functions return arrays of forces along each coordinate axis at each point in time. As explained in :ref:`ref-quickstart`, when observing single or otherwise isolated particles, this is done with *calculate*, and with *calculate_axial* if observing forces acting on a pair of particles. In the latter case, the axes are by default oriented along and perpendicular to the line connecting both traps. In most cases the difference between the inter-particle and inter-trap direction should be rather small, but setting the flag *inter_particle* to *True* as an optional argument switches the axis orientation to the former.

Examples below show both cases. The first three graphs are created by analysing a short sequence of generated data for a single particle, while the latter two are from part of a recorded interaction between two particles.

.. plot:: ../examples/forces_example.py

Use *displacement_calculation* if you need to quickly determine how particles are moving in relation to their traps, for example if the stiffnesses haven't been determined yet. The estimated mean displacement and uncertainties are mostly intended for use with fixed traps.
