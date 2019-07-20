.. _ref-calibration:

Calibration
===========

Let us take a more detailed look into the calibration functionality.

Suppose the potential of an optical trap is elliptical in shape with its long (a) and short (b) axes not, in general, aligned with the coordinate axes of the laboratory sistem. Usually the beam is symmetrical enough for this approximation to make sense. The *calibrate* function then computes the angle :math:`\phi` for which its long axis is rotated with respect to the x-axis and the trap coefficients :math:`k_a, k_b` which characterise the optical trapping force in a parabolic potential, namely :math:`F_i = k_i \Delta r_i`. The new coefficients are defined as

.. math::

    \left [ \matrix{k_a \cr k_b} \right ] = \left [ \matrix{\cos \phi & - \sin \phi \cr \sin \phi & \cos \phi} \right ] \left [ \matrix{k_x \cr k_y} \right ]

and are, together with the rotation angle, the main information about the trap needed for force calculation. Of course, if the displacements from the trap's center are large enough that the parabolic approximation introduces significant error, the shape of the potential can be computed by creating histograms of particle positions. We do not currently provide methods which make use of this data, but since the fitting functions return a list of estimated potential magnitudes at different displacements, all the inputs for adding such are provided.

Coordinate drift is minimised using averaging over time intervals: by default the length is 1 second, but depending on your frequency of measurement you should take care that each average is computed with a large enough number of data points. 

In the following example we can see what the averaging, drift elimination and potential estimation look like in practice - although the plots in this particular example are created on the basis of generated synthetic data, from which trap parameters are computed. First, particle positions are averaged in time; both to better show the movement trend and to subtract it. The effect of this can be seen in the second plot. The final two plots show calculated potential shapes.
    
.. plot:: ../examples/calibration_example.py

