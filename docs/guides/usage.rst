.. _ref-quickstart:

Usage primer
============

Near the focus of a laser beam, particles with refractive indices greater than the surrounding medium are acted upon by an attractive optical gradient force. Optical tweezers exploit this for trapping, and manipulation, of these particles. The purpose of this package is to provide a means of analysing the movement of trapped particles and, consequently, the forces which traps exert on them.

Every important step of the process of data extraction and analysis is included, but only briefly. For more information about the different functionalities, you should consult the respective pages (e.g. :ref:`ref-data-extraction`, :ref:`ref-calibration` etc.). In this quick-start guide, we will assume that experimental data is *not* available. Readers who wish to skip straight to manipulating the data can start with section :ref:`ref-quickstart-calibration`.


.. _ref-generating-sample:

Generating a data sample
------------------------

If measurements are currently unavailable or you want to run tests on a predictable set of data, the *generate* functionality can simulate the behaviour of single particles in static optical traps; to use it, you will need to provide at least the trap strength, given as a tuple of coefficients for the x and y axes. Other optional arguments can change the number of generated points and their spacing in time, decenter the trap and rotate its axes. If the default arguments are unchanged, this works as follows:

.. code-block:: python
    
    times = tweezer.calibration_generate_data.generate_time()
    trajectory = tweezer.calibration_generate_data.generate((kx,ky))
    
The function draws at random from a bivariate distribution and scales the offsets according to the trap coefficients. The x- and y- positions of the particle are returned in a two-column array, while the trap position is of course just the specified *(x,y)* offset, with the default position being centered at (0,0).

If the flag *Drift* is set to *True*, then the particle positions will be convolved with a random sinusoidal function to approximate the effects of random trap position drift due to outside factors.


.. _ref-read-data:

Reading the data
----------------

You can unpack a generated text file with the function *read_file*. Let's say the video we were working with contains two particles:

.. code-block:: python

    import tweezer
    
    time,traps,trajectories = tweezer.plotting.read_file("example_data.dat",2)

The function returns arrays of times, all the trap data and the two particles' trajectories. We save the output and can then analyse it in several ways.


.. _ref-quickstart-calibration:

Calibration
-----------


Even when the exact locations and strengths of the tweezer's optical traps are specified, they do not, in practice, correspond exactly to reality. The laser beam may be slightly offset and not precisely circular; the shape of the trap's potential may not be parabolic even at small displacements or its strength coefficients different than previously computed. This package provides some functionality which reduces the worst issues.

.. note::
    To make the most of the following routines, the particle should only be affected by it's own trap, namely as far removed as possible from any other traps or material it may interact with, including other particles.
    
The *calibrate* function subtracts coordinate drift by averaging positions over time intervals and also rotates the coordinate axes so they are aligned with the short and long axis of the trap's potential shape, which is assumed to be elliptic. It then calculates the trap strength coefficients from the variance of the particle's displacement from trap center.
    
.. code-block:: python

    coefficients,angle,trajectory = tweezer.calibration.calibrate(times,trajectory)
    
Note that it is the particle's, not the trap's trajectory which is altered here, but this should not matter for single-particle examples - but if you are observing interactions between multiple particles or traps, take care when using this function. Rotation of coordinates needs to be manually applied to all other data in order for it to stay correct!

If you are confident that there is no coordinate drift and all you need are mean offsets of the trap center with respect to the specified position, you can use something like the following code, which does not compute the coefficients or rotate data, but it alters the trap positions instead of the particles'.

.. code-block:: python

    offset_x,offset_y = tweezer.offset.single_particle_offset(trap,trajectory)
    trap[:,0] = [x - offset_x for x in trap]
    trap[:,1] = [y - offset_y for y in trap]
    
The offsets of traps at different positions of your workspace can differ - but recording the movements of an unpeturbed trapped particle every time the trap is moved is time-consuming. Assuming that the offset is a linear function of (x,y) coordinates, one identical measurement can be made in each corner of the workspace, from which the offsets at chosen points are interpolated. Let's say the measurements in corners are saved to files named *upper_left.dat*, *upper_right.dat* etc. An example of calculating the offset in arbitrary points (x1,y1) and (y2,y2) is:

.. code-block:: python
    
    points = tweezer.offset.four_corner_offsets("upper_left.dat","upper_right.dat","lower_left.dat","lower_right.dat")
    
    offset_x1,offset_y1 = tweezer.offset.four_corner_calibration(x1, y1, points)
    trap1[:,0] = [x1 - offset_x1 for x in trap1]
    trap1[:,1] = [y1 - offset_y1 for x in trap1]
    
    offset_x2,offset_y2 = tweezer.offset.four_corner_calibration(x2, y2, points)
    trap2[:,0] = [x2 - offset_x2 for x in trap2]
    trap2[:,1] = [y2 - offset_y2 for x in trap2]
    
Analysis of displacements and forces
------------------------------------

With trap coefficients determined and offsets subtracted, forces acting on particles and their displacements can be analysed by functions in *tweezer.force_calc*.

.. code-block:: python
    
    mean_displacements,variances = tweezer.force_calc.displacement_calculation(trajectories, traps)

This returns the mean displacements (in x- and y- directions) of an arbitrary number of particles, as well as their standard deviations. The latter are calculated as if following Poissonian statistics.

If trap coefficients are specified by a (kx,ky) tuple, forces acting on each particle in x- and y-directions can be determined as follows:

.. code-block:: python

    forces,mean_forces = tweezer.force_calc.force_calculation(times,trajectory,trap,coefficients)
    
If the objective is to calculate forces between a *pair* of particles, you can either call the above function twice and manually determine what portion of the force is due to the interaction based on particle positions, or you can use the below method:

.. code-block:: python

    forces,mean_forces,distance = tweezer.force_calc.force_calculation_axis(times, trajectories, traps, coeff_1, coeff_2)
    
    tweezer.plotting.force_plot(times,forces)
    
The forces and means returned are calculated in the axial direction (in the direction away from the other particle), so make sure the sign is correct when calculating the forces' sum. Here, we used *force_plot* at the end to produce a graph of forces with respect to time.


Data format
-----------

When faced with experimental recordings in the form of image sequences or video files, the first course of action is to identify the particles and extract their movement into text files, which can then be further manipulated. This can be done either with TrackPy or another program of your choice. The former can be run through a GUI or directly, which is useful for batch processing.

The tracking packages' final output for each recording is a multi-column .dat file with lines consisting of:

* **[column 1]** time
* **[column 2]** laser power
* **[columns 3-5]** x-, y- coordinate of 1st optical trap and relative trap strength
* same 3 columns repeated for traps 2-4
* **[columns 15-16]** x-,y- coordinates of 1st particle
* same 2 columns repeated for any additional particles

.. note::
    This format only supports up to four optical traps. This was sufficient for our purposes, but in case your measurements require a higher number, you will need to change the read/write functionality.
    

Generating a data sample - moving trap
--------------------------------------

For simulating particles in a moving trap, you can use the slightly more involved (and slower-running) function found in the *data_generation* module. In addition to the trap coefficients, you will need to provide information on the trap movements and viscosity of the sample medium. The latter is due to the fact that here, the virtual particle undergoes stochastic movement between set points in time where its coordinates are recorded. At each (internal) timestep, forces acting on the particle are calculated and its motion updated, although any inter-particle forces are not currently included. This function is called in the following way:

.. code-block:: python

    times,traps,trajectory = tweezer.data_generation.generate(num_points, dt,(kx,ky),trap_frequencies,trap_amplitudes,radius,eta)
    
where trap frequencies and amplitudes are tuples of their values in the x and y directions. By default, the motion along each axis is sinusoidal, although passing *motion_type = 2* as an argument switches to a linear motion.

.. note::
    If simulating both sinusoidal trap motion and random drift, make sure that the amplitude of the former is large enough to dominate over random drift effects.
