==================================
Documentation of library functions
==================================

Docstrings of included functions, roughly grouped by functionality.

Main data generation functions
------------------------------

.. note::
    Do we need to include *brownian.py*?

.. automodule:: tweezer.calibration_generate_data
    :members: draw, generate, generate_time

.. automodule:: tweezer.data_generation
    :members: draw_positions, generate
    
Auxiliary generation functions
------------------------------

.. automodule:: tweezer.calibration_generate_data
    :members: rotate_and_decenter, drift
    
.. automodule:: tweezer.data_generation
    :members: decenter, rotate, drift
    
Calibration
-----------

.. automodule:: tweezer.calibration
    :members:

.. automodule:: tweezer.offset
    :members:
    
Force and displacement calculation
----------------------------------

.. automodule:: tweezer.force_calc
    :members: calculate, calculate_axial, evaluate_force, rot_matrix, rot_matrix_axial, displacement_calculation

Visualisation (plots)
---------------------

.. automodule:: tweezer.plotting
    :members: trajectory_plot, displacement_plot, calibration_plots, potential_plot, potential_polynomial_fit_plot, force_plot, force_plot_radial, force_distance_plot
    
    
Input/Output functionality
==========================

Writing frames
--------------

.. automodule:: tweezer.generate_gaussian_particles
    :members: make_frame, make_frame_array
    
Writing frames - examples
-------------------------

.. automodule:: tweezer.generate_gaussian_particles
    :members: example,two_particles_approaching,two_particles_yield
    
Particle tracking
-----------------

.. automodule:: tweezer.tracking
    :members:
    
File parser
-----------

.. automodule:: tweezer.plotting
    :members: read_file
    
